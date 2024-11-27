// Package multi_csv allows conversion from multi config csv to json.
//
// It is implemented as a go binding to the rust code.
package multi_csv

/*
long long convert_multi_config(
	char *buffer,
	unsigned long buf_len,
	const char *input,
	unsigned long input_len
);
*/
import "C"
import (
	"fmt"
	"runtime"
	"unsafe"
)

// convertMultiConfig wraps convert_multi_config to avoid repeating
// the type-casting in multiple places.
func convertMultiConfig(buf, input []byte) int {
	r := C.convert_multi_config(
		(*C.char)(unsafe.Pointer(unsafe.SliceData(buf))), C.ulong(len(buf)),
		(*C.char)(unsafe.Pointer(unsafe.SliceData(input))), C.ulong(len(input)))
	// buf is going to be kept alive by whatever called this wanting to read
	// its contents after the call, but we do need to make sure that `input`
	// isn't GC'ed in the middle of the call.
	runtime.KeepAlive(input)
	return int(r)
}

// ConvertMultiConfig parses the input as a multi config CSV.
// Returns the JSON representation of the parsed config.
func ConvertMultiConfig(input []byte) ([]byte, error) {
	var buf []byte
	// To avoid unnecessary back and forth with the C API, try to guestimate
	// how much buffer space is likely to be required.
	// Always allocate at least enough to usually get the error message.
	if l := len(input) * 3; l > 256 {
		buf = make([]byte, l)
	} else {
		buf = make([]byte, 256)
	}
	r := convertMultiConfig(buf, input)
	if r == 0 {
		return nil, nil
	}
	if r < 0 {
		l := -r
		if l > len(buf) {
			buf = make([]byte, l)
			r = convertMultiConfig(buf, input)
			l = -r
		}
		return nil, fmt.Errorf(
			"conversion error: %s", buf[:l])
	}
	if r > len(buf) {
		buf = make([]byte, r)
		r = convertMultiConfig(buf, input)
	}
	return buf[:r], nil
}
