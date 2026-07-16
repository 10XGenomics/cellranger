package collector

import (
	"context"
	"encoding/json"
	"os"
)

type envExtractor struct {
	simpleExtractor
	varName string
}

func (e *envExtractor) MarshalJSON() ([]byte, error) {
	b, err := json.Marshal(e.varName)
	if err != nil {
		return nil, err
	}
	return append(append([]byte(`{"env":`), b...), '}'), nil
}

func (e *envExtractor) Start(_ context.Context,
	_ ValueContext, _ *FileContext, result chan<- any) {
	defer close(result)
	_, ok := os.LookupEnv(e.varName)
	result <- ok
}

func MakeEnv(varName string) ValueExtractor {
	return &envExtractor{varName: varName}
}
