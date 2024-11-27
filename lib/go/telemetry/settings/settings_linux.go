package settings

import (
	"fmt"
	"path/filepath"
)

func executableForPid(pid int) (string, error) {
	return filepath.EvalSymlinks(fmt.Sprint("/proc/", pid, "/exe"))
}
