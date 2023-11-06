package fastq

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"syscall"
)

const GZIP_COMPRESSION_LEVEL = 6

func RC_NUCS(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	case 'N':
		return 'N'
	}
	return 0
}

type FastqRecord struct {
	Qname []byte
	Seq   []byte
	Qual  []byte
}

type FastqReader struct {
	source io.ReadCloser
	Reader *bufio.Reader
	Rc     bool
}

type zipFastqReader struct {
	source io.ReadCloser
	cmd    *exec.Cmd
}

// Read up to len(data) bytes.  This is only used as a source for bufio.Reader,
// which knows how to do the repeated reads thing itself - there's no need for
// that here.
func (self *zipFastqReader) Read(data []byte) (int, error) {
	return self.source.Read(data)
}

func (self *zipFastqReader) Close() error {
	err1 := self.source.Close()
	if err := self.cmd.Wait(); err != nil {
		return err
	}
	return err1
}

func newZipFastqReader(path string) (*zipFastqReader, error) {
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
		_ = stdout.Close()
		return nil, err
	}

	return &zipFastqReader{stdout, cmd}, nil
}

func ReverseComplementSeq(seq []byte) []byte {
	newSeq := make([]byte, len(seq))
	for i, nuc := range seq {
		newSeq[len(seq)-i-1] = RC_NUCS(nuc)
	}
	return newSeq
}

func NewFastqReader(path string, rc bool) (*FastqReader, error) {
	var src io.ReadCloser
	var err error
	if strings.HasSuffix(path, ".gz") {
		src, err = newZipFastqReader(path)
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
	if len(qname) == 0 || qname[0] != '@' {
		err := fmt.Errorf(
			"Invalid FASTQ line, expected line with format @<name>: %s",
			qname)
		return err
	}
	record.Qname = qname[1 : len(qname)-1]

	if self.Rc {
		seq, err := self.Reader.ReadSlice('\n')
		if err != nil {
			return err
		}
		record.Seq = ReverseComplementSeq(seq[:len(seq)-1])
	} else {
		seq, err := self.Reader.ReadBytes('\n')
		if err != nil {
			return err
		}
		record.Seq = seq[:len(seq)-1]
	}

	// Read and discard the next line.
	if b, err := self.Reader.ReadSlice('\n'); err != nil && err != bufio.ErrBufferFull {
		return err
	} else if len(b) == 0 || b[0] != '+' {
		return fmt.Errorf(
			"Invalid FASTQ line, expected line with format +, but had prefix %s",
			b)
	} else {
		// discard remainder of line, if any.
		for err == bufio.ErrBufferFull {
			_, err = self.Reader.ReadSlice('\n')
		}
		if err != nil {
			return err
		}
	}

	qual, err := self.Reader.ReadBytes('\n')
	if err != nil {
		return err
	}
	record.Qual = qual[:len(qual)-1]

	return nil
}

func (self *FastqReader) Close() error {
	return self.source.Close()
}

type ZipFastqWriter struct {
	dest   io.WriteCloser
	cmd    *exec.Cmd
	Writer *bufio.Writer
	Path   string
}

func NewZipFastqWriter(path string, mode string) (*ZipFastqWriter, error) {
	self := ZipFastqWriter{
		Path: path,
	}

	var flag int
	switch mode {
	case "w":
		flag = os.O_WRONLY | os.O_CREATE | os.O_TRUNC
	case "a":
		flag = os.O_APPEND | os.O_WRONLY
	default:
		log.Fatal("Invalid open file mode: ", mode)
	}
	dest, err := os.OpenFile(self.Path, flag, 0666)
	if err != nil {
		return nil, err
	}
	// Once the subprocess starts up, we no longer need to have this file open.
	defer dest.Close()

	cmd := exec.Command("gzip", "-c")
	cmd.Stdout = dest
	self.dest, err = cmd.StdinPipe()
	if err != nil {
		return nil, err
	}

	self.Writer = bufio.NewWriter(self.dest)
	self.cmd = cmd
	return &self, cmd.Start()
}

func (self *ZipFastqWriter) write(data []byte) (int, error) {
	// bufio.Writer guarantees it will either write all of the bytes
	// or return an error.
	return self.Writer.Write(data)
}

func (self *ZipFastqWriter) writeByte(c byte) error {
	return self.Writer.WriteByte(c)
}

func (self *ZipFastqWriter) writeString(data string) (int, error) {
	return self.Writer.WriteString(data)
}

func (self *ZipFastqWriter) WriteRecord(record *FastqRecord) (int, error) {
	if err := self.writeByte('@'); err != nil {
		return 0, err
	}
	if _, err := self.write(record.Qname); err != nil {
		return 0, err
	}
	if err := self.writeByte('\n'); err != nil {
		return 0, err
	}
	if _, err := self.write(record.Seq); err != nil {
		return 0, err
	}
	if _, err := self.writeString("\n+\n"); err != nil {
		return 0, err
	}
	if n, err := self.write(record.Qual); err != nil {
		return n, err
	} else if err := self.writeByte('\n'); err != nil {
		return n, err
	} else {
		return n + 1, err
	}
}

func (self *ZipFastqWriter) Close() error {
	err := self.Writer.Flush()
	if cerr := self.dest.Close(); cerr != nil && err == nil {
		err = cerr
	}
	if cerr := self.cmd.Wait(); cerr != nil && err == nil {
		err = cerr
	}
	return err
}

type ZipFastqWriterCache struct {
	HaveOpened map[string]struct{}
	OpenFiles  OrderedMap
	MaxFiles   int
}

func NewZipFastqWriterCache() (*ZipFastqWriterCache, error) {
	self := &ZipFastqWriterCache{}
	return self, self.Init()
}

func (self *ZipFastqWriterCache) Init() error {
	if err := self.configMaxFiles(); err != nil {
		return err
	}
	self.HaveOpened = make(map[string]struct{}, 128)
	self.OpenFiles.Init(128)
	return nil
}

func (self *ZipFastqWriterCache) configMaxFiles() error {
	var rLimit syscall.Rlimit
	if err := syscall.Getrlimit(syscall.RLIMIT_NOFILE, &rLimit); err != nil {
		return err
	}
	// Leave some headroom for other stuff, but also we have 2 handles for each
	// file: one for the on-disk file, and the other for the pipe to the
	// subprocess.
	self.MaxFiles = (int(rLimit.Cur) - 100) / 2
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
		self.HaveOpened[closeFilename] = struct{}{}
		if err := closeZipFastqWriter.(*ZipFastqWriter).Close(); err != nil {
			return nil, err
		}
	}

	zipFastqWriter, err := NewZipFastqWriter(filename, mode)
	if err != nil {
		return nil, err
	}

	self.OpenFiles.Set(filename, zipFastqWriter)
	return zipFastqWriter, nil
}

func (self *ZipFastqWriterCache) Close() error {
	var lastErr error
	for _, value, ok := self.OpenFiles.Pop(); ok; _, value, ok = self.OpenFiles.Pop() {
		if err := value.(*ZipFastqWriter).Close(); err != nil {
			lastErr = err
		}
	}
	return lastErr
}

type DemultFastqWriter struct {
	SampleIndexes      map[string]struct{}
	Cache              ZipFastqWriterCache
	Path               string
	suffix             string
	DemultReadType     string
	OtherIndexReadType string
	prefix             []byte
	Lane               int
	Chunk              int
}

func NewDemultFastqWriter(path, demultReadType, otherIndexReadType string, lane, chunk int,
	commonSampleIndexes []string) (*DemultFastqWriter, error) {
	self := &DemultFastqWriter{
		Path:               filepath.Clean(path),
		DemultReadType:     demultReadType,
		OtherIndexReadType: otherIndexReadType,
		Lane:               lane,
		Chunk:              chunk,
		SampleIndexes:      make(map[string]struct{}, len(commonSampleIndexes)),
	}
	self.suffix = fmt.Sprintf("_lane-%03d-chunk-%03d.fastq.gz", self.Lane, self.Chunk)
	// Reusable buffer for constructing paths.
	self.prefix = append(append(
		make([]byte, 0, len(self.Path)+len("/read-RA_si-")+16+len(self.suffix)),
		self.Path...), "/read-"...)

	// Initialize file handle cache
	if err := self.Cache.Init(); err != nil {
		return nil, err
	}

	// Initialize sample indices
	for _, commonSampleIndex := range commonSampleIndexes {
		self.SampleIndexes[commonSampleIndex] = struct{}{}
	}

	return self, os.MkdirAll(self.Path, 0755)
}

// Generate the path for the fastq corresponding to the given sequence
// string and read type.
//
// Not thread safe!  Reuses a buffer, transiently, for building the string.
func (self *DemultFastqWriter) makeFastqPath(seq []byte, readType string) string {
	// We run this an awful lot.  Use a reusable buffer to avoid excessive
	// allocations.
	return string(append(append(append(append(
		self.prefix, readType...), "_si-"...), seq...), self.suffix...))
}

// Write records to the appropriate fastq files.
//
// Not thread safe!  Calling from multiple threads can lead to file corruption.
func (self *DemultFastqWriter) WriteRecords(siRecord *FastqRecord,
	read1Record *FastqRecord, read2Record *FastqRecord,
	otherIndexRecord *FastqRecord, siName []byte) error {
	seq := INVALID_SAMPLE_INDEX
	// Override with siName, or sampleIndex seq if valid
	if len(siName) > 0 {
		seq = siName
	} else if _, ok := self.SampleIndexes[string(siRecord.Seq)]; ok {
		seq = siRecord.Seq
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

func (self *DemultFastqWriter) Close() error {
	return self.Cache.Close()
}
