package converter

import (
	"encoding/json"
	"fmt"

	"github.com/10XDev/cellranger/lib/go/multi_csv"
)

type Converter interface {
	json.Marshaler
	Convert(b []byte) ([]byte, error)
}

type multicsv struct{}

func (multicsv) Convert(b []byte) ([]byte, error) {
	return multi_csv.ConvertMultiConfig(b)
}

// MarshalJSON implements Converter.
func (multicsv) MarshalJSON() ([]byte, error) {
	return []byte(`{"converter":"multicsv"}`), nil
}

func Make(name string) (Converter, error) {
	if name == "" {
		return nil, nil
	}
	switch name {
	case "multicsv":
		return multicsv{}, nil
	}
	return nil, fmt.Errorf("converter %q not implemented", name)
}
