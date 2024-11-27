package collector

// #include <gnu/libc-version.h>
import "C"
import "context"

type glibcExtractor struct{ simpleExtractor }

func (glibcExtractor) MarshalJSON() ([]byte, error) {
	return []byte(`{"special":"glibc"}`), nil
}

// Start implements ValueExtractor.
func (glibcExtractor) Start(_ context.Context, _ ValueContext, _ *FileContext, result chan<- any) {
	defer close(result)
	result <- C.GoString(C.gnu_get_libc_version())
}
