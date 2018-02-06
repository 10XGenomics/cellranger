package fastq

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/youtube/vitess/go/cgzip"
	"io"
	"log"
	"os"
	"os/exec"
	"path"
	"strings"
	"syscall"
)

const GZIP_COMPRESSION_LEVEL = 6

var RC_NUCS map[byte]byte

type FastqRecord struct {
	Qname []byte
	Seq   []byte
	Qual  []byte
}

type FastqReader struct {
	Source io.ReadCloser
	Reader *bufio.Reader
	Rc     bool
}

type ZipFastqReader struct {
	Source io.ReadCloser
	Cmd    *exec.Cmd
}

func (self *ZipFastqReader) Read(data []byte) (int, error) {
	want := len(data)
	var err error
	var offset int
	var read_len int
	for offset = 0; offset < want && err == nil; read_len, err = self.Source.Read(data[offset:]) {
		offset += read_len
	}

	return offset, err
}

func (self *ZipFastqReader) Close() error {
	self.Source.Close()
	return self.Cmd.Wait()
}

func NewZipFastqReader(path string) (*ZipFastqReader, error) {
	if _, err := os.Stat(path); os.IsNotExist(err) {
		panic("no such file or directory: " + path)
	}
	cmd := exec.Command("gunzip", "-c", path)
	stdout, err := cmd.StdoutPipe()

	if err != nil {
		return nil, err
	}
	err = cmd.Start()

	if err != nil {
		return nil, err
	}

	return &ZipFastqReader{stdout, cmd}, nil
}

func ReverseComplementSeq(seq []byte) []byte {
	newSeq := make([]byte, len(seq))
	for i, nuc := range seq {
		newSeq[len(seq)-i-1] = RC_NUCS[nuc]
	}
	return newSeq
}

func NewFastqReader(path string, rc bool) (*FastqReader, error) {
	RC_NUCS = map[byte]byte{
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'N': 'N',
	}

	var src io.ReadCloser
	var err error
	if strings.HasSuffix(path, ".gz") {
		src, err = NewZipFastqReader(path)
	} else {
		src, err = os.Open(path)
	}
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(src)
	return &FastqReader{src, reader, rc}, nil
}

func (self *FastqReader) ReadRecord(record *FastqRecord) error {
	qname, err := self.Reader.ReadBytes('\n')
	if err != nil {
		return err
	}
	if qname[0] != '@' {
		err := errors.New(fmt.Sprintf("Invalid FASTQ line, expected line with format @<name>: %s", qname))
		return err
	}
	record.Qname = qname[1 : len(qname)-1]

	seq, err := self.Reader.ReadBytes('\n')
	if err != nil {
		return err
	}
	record.Seq = seq[0 : len(seq)-1]
	if self.Rc {
		record.Seq = ReverseComplementSeq(record.Seq)
	}

	line, err := self.Reader.ReadBytes('\n')
	if err != nil {
		return err
	}
	if line[0] != '+' {
		err := errors.New(fmt.Sprintf("Invalid FASTQ line, expected line with format +: %s", line))
		return err
	}

	qual, err := self.Reader.ReadBytes('\n')
	if err != nil {
		return err
	}
	record.Qual = qual[0 : len(qual)-1]

	return nil
}

func (self *FastqReader) Close() error {
	return self.Source.Close()
}

type ZipFastqWriter struct {
	Path             string
	Source           io.WriteCloser
	CompressedSource *cgzip.Writer
	Writer           *bufio.Writer
}

func NewZipFastqWriter(path string, mode string) (*ZipFastqWriter, error) {
	self := &ZipFastqWriter{}
	self.Path = path

	var err error
	if mode == "w" {
		self.Source, err = os.Create(self.Path)
	} else if mode == "a" {
		self.Source, err = os.OpenFile(self.Path, os.O_APPEND|os.O_WRONLY, 0755)
	} else {
		log.Fatal("Invalid open file mode: ", mode)
	}
	if err != nil {
		return nil, err
	}

	self.CompressedSource, err = cgzip.NewWriterLevel(self.Source, GZIP_COMPRESSION_LEVEL)
	if err != nil {
		return nil, err
	}

	self.Writer = bufio.NewWriter(self.CompressedSource)
	return self, nil
}

func (self *ZipFastqWriter) write(data []byte) (int, error) {
	want := len(data)
	var err error
	var offset int
	var write_len int
	for offset = 0; offset < want && err == nil; write_len, err = self.Writer.Write(data[offset:]) {
		offset += write_len
	}

	return offset, err
}

func (self *ZipFastqWriter) WriteRecord(record *FastqRecord) (int, error) {
	self.write([]byte{'@'})
	self.write(append(record.Qname, '\n'))
	self.write(append(record.Seq, '\n'))
	self.write([]byte{'+', '\n'})
	return self.write(append(record.Qual, '\n'))
}

func (self *ZipFastqWriter) Close() error {
	self.Writer.Flush()
	self.CompressedSource.Close()
	return self.Source.Close()
}

type ZipFastqWriterCache struct {
	MaxFiles   int
	HaveOpened map[string]bool
	OpenFiles  *OrderedMap
}

func NewZipFastqWriterCache() (*ZipFastqWriterCache, error) {
	self := &ZipFastqWriterCache{}
	if err := self.configMaxFiles(); err != nil {
		return nil, err
	}
	self.HaveOpened = map[string]bool{}
	self.OpenFiles = NewOrderedMap()
	return self, nil
}

func (self *ZipFastqWriterCache) configMaxFiles() error {
	var rLimit syscall.Rlimit
	if err := syscall.Getrlimit(syscall.RLIMIT_NOFILE, &rLimit); err != nil {
		return err
	}
	self.MaxFiles = int(rLimit.Cur) - 100
	if self.MaxFiles < 0 {
		log.Fatal("Soft file handle limit is too low: ", rLimit.Cur)
	}
	return nil
}

func (self *ZipFastqWriterCache) Get(filename string) (*ZipFastqWriter, error) {
	if self.OpenFiles.HasKey(filename) {
		zipFastqWriter, _ := self.OpenFiles.Get(filename)
		self.OpenFiles.Set(filename, zipFastqWriter)
		return zipFastqWriter.(*ZipFastqWriter), nil
	}

	var mode string
	if _, ok := self.HaveOpened[filename]; ok {
		mode = "a"
	} else {
		mode = "w"
	}

	if self.OpenFiles.Length() == self.MaxFiles-1 {
		closeFilename, closeZipFastqWriter, _ := self.OpenFiles.Pop()
		self.HaveOpened[closeFilename] = true
		closeZipFastqWriter.(*ZipFastqWriter).Close()
	}

	zipFastqWriter, err := NewZipFastqWriter(filename, mode)
	if err != nil {
		return nil, err
	}

	self.OpenFiles.Set(filename, zipFastqWriter)
	return zipFastqWriter, nil
}

func (self *ZipFastqWriterCache) Close() {
	for _, value := range self.OpenFiles.Values() {
		value.(*ZipFastqWriter).Close()
	}
}

type DemultFastqWriter struct {
	SampleIndexes      map[string]bool
	Cache              *ZipFastqWriterCache
	Path               string
	DemultReadType     string
	OtherIndexReadType string
	Lane               int
	Chunk              int
}

func NewDemultFastqWriter(path string, demultReadType string, otherIndexReadType string, lane int, chunk int,
	commonSampleIndexes []string) (*DemultFastqWriter, error) {
	self := &DemultFastqWriter{}
	self.Path = path
	self.DemultReadType = demultReadType
	self.OtherIndexReadType = otherIndexReadType
	self.Lane = lane
	self.Chunk = chunk

	// Initialize file handle cache
	var err error
	self.Cache, err = NewZipFastqWriterCache()
	if err != nil {
		return nil, err
	}

	// Create directory
	os.MkdirAll(self.Path, 0755)

	// Initialize sample indices
	self.SampleIndexes = map[string]bool{}
	for _, commonSampleIndex := range commonSampleIndexes {
		self.SampleIndexes[commonSampleIndex] = true
	}

	return self, nil
}

func (self *DemultFastqWriter) makeFastqPath(seq string, readType string) string {
	return path.Join(self.Path, fmt.Sprintf("read-%s_si-%s_lane-%03d-chunk-%03d.fastq.gz", readType, seq, self.Lane, self.Chunk))
}

func (self *DemultFastqWriter) WriteRecords(siRecord *FastqRecord, read1Record *FastqRecord, read2Record *FastqRecord, otherIndexRecord *FastqRecord) error {
	seq := string(siRecord.Seq)
	if _, ok := self.SampleIndexes[seq]; !ok {
		seq = INVALID_SAMPLE_INDEX
	}

	// Write sample index
	siFastqPath := self.makeFastqPath(seq, self.DemultReadType)
	f, err := self.Cache.Get(siFastqPath)
	if err != nil {
		return err
	}
	if _, err := f.WriteRecord(siRecord); err != nil {
		return err
	}

	// Write read 1
	readFastqPath := self.makeFastqPath(seq, INTERLEAVED_READ_TYPE)
	f, err = self.Cache.Get(readFastqPath)
	if err != nil {
		return err
	}
	if _, err := f.WriteRecord(read1Record); err != nil {
		return err
	}

	// Write read 2
	if len(read2Record.Seq) > 0 {
		if _, err := f.WriteRecord(read2Record); err != nil {
			return err
		}
	}

	// Write other index
	if len(otherIndexRecord.Seq) > 0 {
		otherIndexFastqPath := self.makeFastqPath(seq, self.OtherIndexReadType)
		f, err = self.Cache.Get(otherIndexFastqPath)
		if err != nil {
			return err
		}
		if _, err := f.WriteRecord(otherIndexRecord); err != nil {
			return err
		}
	}

	return nil
}

func (self *DemultFastqWriter) Close() {
	self.Cache.Close()
}
