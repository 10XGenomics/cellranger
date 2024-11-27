package collector

import (
	"context"
	"encoding/json"
	"os"
	"strings"
)

type flagExtractor struct {
	simpleExtractor
	flag string
	mro  bool
}

func (f *flagExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(f.flag)
	if err != nil {
		return nil, err
	}
	if f.mro {
		return append(append([]byte(`{"mro_flag":`), b...), '}'), nil
	} else {
		return append(append([]byte(`{"flag":`), b...), '}'), nil
	}
}

func (f *flagExtractor) Start(_ context.Context,
	vc ValueContext, _ *FileContext, result chan<- any) {
	defer close(result)
	if f.mro {
		for _, flag := range vc.MroFlags() {
			if before, after, found := strings.Cut(flag, "="); found && before == f.flag {
				result <- after
				return
			}
		}
		for _, flag := range strings.Fields(os.Getenv("MROFLAGS")) {
			if before, after, found := strings.Cut(flag, "="); found && before == f.flag {
				result <- after
				return
			}
		}
		return
	}
	for _, flag := range vc.CommandLine() {
		if before, after, found := strings.Cut(flag, "="); found && before == f.flag {
			result <- after
			return
		}
	}
}

func MakeFlag(flag string, mro bool) ValueExtractor {
	return &flagExtractor{
		flag: flag,
		mro:  mro,
	}
}
