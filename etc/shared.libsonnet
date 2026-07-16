// Shared telemetry constants used across cellranger, spaceranger, and cellranger-arc.
// Do not edit manually; regenerated from the golden telemetry JSON outputs.
local command_group = { when: 'command' };
local subcommand_group(cmd, when) = { subcommands: [cmd], when: when };
{
  // Common constructor for metrics where name == json_path and values.float is true.
  path_metric_float(name, path, groups): {
    file_content: {
      json_path: [
        name,
      ],
      path: [
        path,
      ],
    },
    groups: groups,
    name: name,
    values: {
      float: true,
    },
  },

  // Groups identical across all products.
  groups: {
    aggr_failure: subcommand_group('aggr', 'failure'),
    count_failure: subcommand_group('count', 'failure'),
    hardware_stats: command_group,
    mkfastq_failure: subcommand_group('mkfastq', 'failure'),
    mkfastq_success: subcommand_group('mkfastq', 'success'),
    mkref_failure: subcommand_group('mkref', 'failure'),
    mkref_success: subcommand_group('mkref', 'success'),
    os_stats: command_group,
    reanalyze_failure: subcommand_group('reanalyze', 'failure'),
    reanalyze_success: subcommand_group('reanalyze', 'success'),
  },

  // Error messages
  error_matches:
    // Cell Ranger TXRNGR error codes
    [
      { exp: 'TXRNGR%d' % code }
      for code in std.range(10001, 10128)
    ] +
    [
      // Martian errors
      { key: 'Martian: Job failed in stage code', exp: 'Job failed in stage code' },
      { key: 'Martian: not queued or running', exp: 'not queued or running' },
      // POSIX errors
      { key: 'EPERM', exp: '(?i:Operation not permitted)' },  // errno=1
      { key: 'ENOENT', exp: '(?i:No such file or directory)' },  // errno=2
      { key: 'ESRCH', exp: '(?i:No such process)' },  // errno=3
      { key: 'EINTR', exp: '(?i:Interrupted system call)' },  // errno=4
      { key: 'EIO', exp: '(?i:Input/output error)' },  // errno=5
      { key: 'ENXIO', exp: '(?i:No such device or address)' },  // errno=6
      { key: 'E2BIG', exp: '(?i:Argument list too long)' },  // errno=7
      { key: 'ENOEXEC', exp: '(?i:Exec format error)' },  // errno=8
      { key: 'EBADF', exp: '(?i:Bad file descriptor)' },  // errno=9
      { key: 'ECHILD', exp: '(?i:No child processes)' },  // errno=10
      { key: 'EAGAIN', exp: '(?i:Resource temporarily unavailable)' },  // errno=11
      { key: 'ENOMEM', exp: '(?i:Cannot allocate memory)' },  // errno=12
      { key: 'EACCES', exp: '(?i:Permission denied)' },  // errno=13
      { key: 'EFAULT', exp: '(?i:Bad address)' },  // errno=14
      { key: 'ENOTBLK', exp: '(?i:Block device required)' },  // errno=15
      { key: 'EBUSY', exp: '(?i:Device or resource busy)' },  // errno=16
      { key: 'EEXIST', exp: '(?i:File exists)' },  // errno=17
      { key: 'EXDEV', exp: '(?i:Invalid cross-device link)' },  // errno=18
      { key: 'ENODEV', exp: '(?i:No such device)' },  // errno=19
      { key: 'ENOTDIR', exp: '(?i:Not a directory)' },  // errno=20
      { key: 'EISDIR', exp: '(?i:Is a directory)' },  // errno=21
      { key: 'EINVAL', exp: '(?i:Invalid argument)' },  // errno=22
      { key: 'ENFILE', exp: '(?i:Too many open files in system)' },  // errno=23
      { key: 'EMFILE', exp: '(?i:Too many open files)' },  // errno=24
      { key: 'ENOTTY', exp: '(?i:Inappropriate ioctl for device)' },  // errno=25
      { key: 'ETXTBSY', exp: '(?i:Text file busy)' },  // errno=26
      { key: 'EFBIG', exp: '(?i:File too large)' },  // errno=27
      { key: 'ENOSPC', exp: '(?i:No space left on device)' },  // errno=28
      { key: 'ESPIPE', exp: '(?i:Illegal seek)' },  // errno=29
      { key: 'EROFS', exp: '(?i:Read-only file system)' },  // errno=30
      { key: 'EMLINK', exp: '(?i:Too many links)' },  // errno=31
      { key: 'EPIPE', exp: '(?i:Broken pipe)' },  // errno=32
      { key: 'EDOM', exp: '(?i:Numerical argument out of domain)' },  // errno=33
      { key: 'ERANGE', exp: '(?i:Numerical result out of range)' },  // errno=34
      { key: 'EDEADLK', exp: '(?i:Resource deadlock avoided)' },  // errno=35
      { key: 'ENAMETOOLONG', exp: '(?i:File name too long)' },  // errno=36
      { key: 'ENOLCK', exp: '(?i:No locks available)' },  // errno=37
      { key: 'ENOSYS', exp: '(?i:Function not implemented)' },  // errno=38
      { key: 'ENOTEMPTY', exp: '(?i:Directory not empty)' },  // errno=39
      { key: 'ELOOP', exp: '(?i:Too many levels of symbolic links)' },  // errno=40
      { key: 'ENOMSG', exp: '(?i:No message of desired type)' },  // errno=42
      { key: 'EIDRM', exp: '(?i:Identifier removed)' },  // errno=43
      { key: 'ECHRNG', exp: '(?i:Channel number out of range)' },  // errno=44
      { key: 'EL2NSYNC', exp: '(?i:Level 2 not synchronized)' },  // errno=45
      { key: 'EL3HLT', exp: '(?i:Level 3 halted)' },  // errno=46
      { key: 'EL3RST', exp: '(?i:Level 3 reset)' },  // errno=47
      { key: 'ELNRNG', exp: '(?i:Link number out of range)' },  // errno=48
      { key: 'EUNATCH', exp: '(?i:Protocol driver not attached)' },  // errno=49
      { key: 'ENOCSI', exp: '(?i:No CSI structure available)' },  // errno=50
      { key: 'EL2HLT', exp: '(?i:Level 2 halted)' },  // errno=51
      { key: 'EBADE', exp: '(?i:Invalid exchange)' },  // errno=52
      { key: 'EBADR', exp: '(?i:Invalid request descriptor)' },  // errno=53
      { key: 'EXFULL', exp: '(?i:Exchange full)' },  // errno=54
      { key: 'ENOANO', exp: '(?i:No anode)' },  // errno=55
      { key: 'EBADRQC', exp: '(?i:Invalid request code)' },  // errno=56
      { key: 'EBADSLT', exp: '(?i:Invalid slot)' },  // errno=57
      { key: 'EBFONT', exp: '(?i:Bad font file format)' },  // errno=59
      { key: 'ENOSTR', exp: '(?i:Device not a stream)' },  // errno=60
      { key: 'ENODATA', exp: '(?i:No data available)' },  // errno=61
      { key: 'ETIME', exp: '(?i:Timer expired)' },  // errno=62
      { key: 'ENOSR', exp: '(?i:Out of streams resources)' },  // errno=63
      { key: 'ENONET', exp: '(?i:Machine is not on the network)' },  // errno=64
      { key: 'ENOPKG', exp: '(?i:Package not installed)' },  // errno=65
      { key: 'EREMOTE', exp: '(?i:Object is remote)' },  // errno=66
      { key: 'ENOLINK', exp: '(?i:Link has been severed)' },  // errno=67
      { key: 'EADV', exp: '(?i:Advertise error)' },  // errno=68
      { key: 'ESRMNT', exp: '(?i:Srmount error)' },  // errno=69
      { key: 'ECOMM', exp: '(?i:Communication error on send)' },  // errno=70
      { key: 'EPROTO', exp: '(?i:Protocol error)' },  // errno=71
      { key: 'EMULTIHOP', exp: '(?i:Multihop attempted)' },  // errno=72
      { key: 'EDOTDOT', exp: '(?i:RFS specific error)' },  // errno=73
      { key: 'EBADMSG', exp: '(?i:Bad message)' },  // errno=74
      { key: 'EOVERFLOW', exp: '(?i:Value too large for defined data type)' },  // errno=75
      { key: 'ENOTUNIQ', exp: '(?i:Name not unique on network)' },  // errno=76
      { key: 'EBADFD', exp: '(?i:File descriptor in bad state)' },  // errno=77
      { key: 'EREMCHG', exp: '(?i:Remote address changed)' },  // errno=78
      { key: 'ELIBACC', exp: '(?i:Can not access a needed shared library)' },  // errno=79
      { key: 'ELIBBAD', exp: '(?i:Accessing a corrupted shared library)' },  // errno=80
      { key: 'ELIBSCN', exp: '(?i:\\.lib section in a\\.out corrupted)' },  // errno=81
      { key: 'ELIBMAX', exp: '(?i:Attempting to link in too many shared libraries)' },  // errno=82
      { key: 'ELIBEXEC', exp: '(?i:Cannot exec a shared library directly)' },  // errno=83
      { key: 'EILSEQ', exp: '(?i:Invalid or incomplete multibyte or wide character)' },  // errno=84
      { key: 'ERESTART', exp: '(?i:Interrupted system call should be restarted)' },  // errno=85
      { key: 'ESTRPIPE', exp: '(?i:Streams pipe error)' },  // errno=86
      { key: 'EUSERS', exp: '(?i:Too many users)' },  // errno=87
      { key: 'ENOTSOCK', exp: '(?i:Socket operation on non-socket)' },  // errno=88
      { key: 'EDESTADDRREQ', exp: '(?i:Destination address required)' },  // errno=89
      { key: 'EMSGSIZE', exp: '(?i:Message too long)' },  // errno=90
      { key: 'EPROTOTYPE', exp: '(?i:Protocol wrong type for socket)' },  // errno=91
      { key: 'ENOPROTOOPT', exp: '(?i:Protocol not available)' },  // errno=92
      { key: 'EPROTONOSUPPORT', exp: '(?i:Protocol not supported)' },  // errno=93
      { key: 'ESOCKTNOSUPPORT', exp: '(?i:Socket type not supported)' },  // errno=94
      { key: 'EOPNOTSUPP', exp: '(?i:Operation not supported)' },  // errno=95
      { key: 'EPFNOSUPPORT', exp: '(?i:Protocol family not supported)' },  // errno=96
      { key: 'EAFNOSUPPORT', exp: '(?i:Address family not supported by protocol)' },  // errno=97
      { key: 'EADDRINUSE', exp: '(?i:Address already in use)' },  // errno=98
      { key: 'EADDRNOTAVAIL', exp: '(?i:Cannot assign requested address)' },  // errno=99
      { key: 'ENETDOWN', exp: '(?i:Network is down)' },  // errno=100
      { key: 'ENETUNREACH', exp: '(?i:Network is unreachable)' },  // errno=101
      { key: 'ENETRESET', exp: '(?i:Network dropped connection on reset)' },  // errno=102
      { key: 'ECONNABORTED', exp: '(?i:Software caused connection abort)' },  // errno=103
      { key: 'ECONNRESET', exp: '(?i:Connection reset by peer)' },  // errno=104
      { key: 'ENOBUFS', exp: '(?i:No buffer space available)' },  // errno=105
      { key: 'EISCONN', exp: '(?i:Transport endpoint is already connected)' },  // errno=106
      { key: 'ENOTCONN', exp: '(?i:Transport endpoint is not connected)' },  // errno=107
      { key: 'ESHUTDOWN', exp: '(?i:Cannot send after transport endpoint shutdown)' },  // errno=108
      { key: 'ETOOMANYREFS', exp: '(?i:Too many references: cannot splice)' },  // errno=109
      { key: 'ETIMEDOUT', exp: '(?i:Connection timed out)' },  // errno=110
      { key: 'ECONNREFUSED', exp: '(?i:Connection refused)' },  // errno=111
      { key: 'EHOSTDOWN', exp: '(?i:Host is down)' },  // errno=112
      { key: 'EHOSTUNREACH', exp: '(?i:No route to host)' },  // errno=113
      { key: 'EALREADY', exp: '(?i:Operation already in progress)' },  // errno=114
      { key: 'EINPROGRESS', exp: '(?i:Operation now in progress)' },  // errno=115
      { key: 'ESTALE', exp: '(?i:Stale file handle)' },  // errno=116
      { key: 'EUCLEAN', exp: '(?i:Structure needs cleaning)' },  // errno=117
      { key: 'ENOTNAM', exp: '(?i:Not a Xenix named type file)' },  // errno=118
      { key: 'ENAVAIL', exp: '(?i:No Xenix semaphores available)' },  // errno=119
      { key: 'EISNAM', exp: '(?i:Is a named type file)' },  // errno=120
      { key: 'EREMOTEIO', exp: '(?i:Remote I/O error)' },  // errno=121
      { key: 'EDQUOT', exp: '(?i:Disk quota exceeded)' },  // errno=122
      { key: 'ENOMEDIUM', exp: '(?i:No medium found)' },  // errno=123
      { key: 'EMEDIUMTYPE', exp: '(?i:Wrong medium type)' },  // errno=124
      { key: 'ECANCELED', exp: '(?i:Operation canceled)' },  // errno=125
      { key: 'ENOKEY', exp: '(?i:Required key not available)' },  // errno=126
      { key: 'EKEYEXPIRED', exp: '(?i:Key has expired)' },  // errno=127
      { key: 'EKEYREVOKED', exp: '(?i:Key has been revoked)' },  // errno=128
      { key: 'EKEYREJECTED', exp: '(?i:Key was rejected by service)' },  // errno=129
      { key: 'EOWNERDEAD', exp: '(?i:Owner died)' },  // errno=130
      { key: 'ENOTRECOVERABLE', exp: '(?i:State not recoverable)' },  // errno=131
      { key: 'ERFKILL', exp: '(?i:Operation not possible due to RF-kill)' },  // errno=132
      { key: 'EHWPOISON', exp: '(?i:Memory page has hardware error)' },  // errno=133
      { key: 'ENOTSUP', exp: '(?i:Not supported parameter or option)' },  // errno=134
      // Python exceptions
      self.prefix_re_matcher('AssertionError'),
      self.prefix_re_matcher('AttributeError'),
      self.prefix_re_matcher('EOFError'),
      self.prefix_re_matcher('FloatingPointError'),
      self.prefix_re_matcher('GeneratorExit'),
      self.prefix_re_matcher('ImportError'),
      self.prefix_re_matcher('ModuleNotFoundError'),
      self.prefix_re_matcher('IndexError'),
      self.prefix_re_matcher('KeyError'),
      self.prefix_re_matcher('KeyboardInterrupt'),
      self.prefix_re_matcher('MemoryError'),
      self.prefix_re_matcher('NameError'),
      self.prefix_re_matcher('NotImplementedError'),
      self.prefix_re_matcher('OSError'),
      self.prefix_re_matcher('OverflowError'),
      self.prefix_re_matcher('PythonFinalizationError'),
      self.prefix_re_matcher('RecursionError'),
      self.prefix_re_matcher('ReferenceError'),
      self.prefix_re_matcher('RuntimeError'),
      self.prefix_re_matcher('StopIteration'),
      self.prefix_re_matcher('StopAsyncIteration'),
      self.prefix_re_matcher('SyntaxError'),
      self.prefix_re_matcher('IndentationError'),
      self.prefix_re_matcher('TabError'),
      self.prefix_re_matcher('SystemError'),
      self.prefix_re_matcher('SystemExit'),
      self.prefix_re_matcher('TypeError'),
      self.prefix_re_matcher('UnboundLocalError'),
      self.prefix_re_matcher('UnicodeError'),
      self.prefix_re_matcher('UnicodeEncodeError'),
      self.prefix_re_matcher('UnicodeDecodeError'),
      self.prefix_re_matcher('UnicodeTranslateError'),
      self.prefix_re_matcher('ValueError'),
      self.prefix_re_matcher('ZeroDivisionError'),
      // OS Exceptions (OSError subclasses)
      self.prefix_re_matcher('BlockingIOError'),
      self.prefix_re_matcher('ChildProcessError'),
      self.prefix_re_matcher('ConnectionError'),
      self.prefix_re_matcher('BrokenPipeError'),
      self.prefix_re_matcher('ConnectionAbortedError'),
      self.prefix_re_matcher('ConnectionRefusedError'),
      self.prefix_re_matcher('ConnectionResetError'),
      self.prefix_re_matcher('FileExistsError'),
      self.prefix_re_matcher('FileNotFoundError'),
      self.prefix_re_matcher('InterruptedError'),
      self.prefix_re_matcher('IsADirectoryError'),
      self.prefix_re_matcher('NotADirectoryError'),
      self.prefix_re_matcher('PermissionError'),
      self.prefix_re_matcher('ProcessLookupError'),
      self.prefix_re_matcher('TimeoutError'),
      // Cell Ranger Python exceptions
      self.prefix_re_matcher('AggrDefError'),
      self.prefix_re_matcher('AmbiguousValueError'),
      self.prefix_re_matcher('AnnotationError'),
      self.prefix_re_matcher('ArcReferenceError'),
      self.prefix_re_matcher('AtacGexReferenceError'),
      self.prefix_re_matcher('BarcodeAssignmentException'),
      self.prefix_re_matcher('BaseReferenceError'),
      self.prefix_re_matcher('BedFormatError'),
      self.prefix_re_matcher('CellCallerException'),
      self.prefix_re_matcher('CellEstimateCantConvergeException'),
      self.prefix_re_matcher('CelltypeDefException'),
      self.prefix_re_matcher('CloudAPIError'),
      self.prefix_re_matcher('CloudCellAnnotationException'),
      self.prefix_re_matcher('CometFileNotFoundError'),
      self.prefix_re_matcher('ConfigParseError'),
      self.prefix_re_matcher('CopyFastqsError'),
      self.prefix_re_matcher('CRCalledProcessError'),
      self.prefix_re_matcher('CSVEmptyException'),
      self.prefix_re_matcher('CSVParseException'),
      self.prefix_re_matcher('DualIndexFlowcellException'),
      self.prefix_re_matcher('DuplicateContigNameException'),
      self.prefix_re_matcher('EmptyFileError'),
      self.prefix_re_matcher('FastqError'),
      self.prefix_re_matcher('FastqParseError'),
      self.prefix_re_matcher('FeatureDefException'),
      self.prefix_re_matcher('FeatureLinkageError'),
      self.prefix_re_matcher('ForceCellsError'),
      self.prefix_re_matcher('FragmentException'),
      self.prefix_re_matcher('GexReferenceError'),
      self.prefix_re_matcher('GtfParseError'),
      self.prefix_re_matcher('HeaderException'),
      self.prefix_re_matcher('IDNotFoundError'),
      self.prefix_re_matcher('IndexAmbiguityException'),
      self.prefix_re_matcher('InvalidCellBarcode'),
      self.prefix_re_matcher('InvalidLibraryError'),
      self.prefix_re_matcher('JibesError'),
      self.prefix_re_matcher('MatrixRankTooSmallException'),
      self.prefix_re_matcher('MergeError'),
      self.prefix_re_matcher('MetricsError'),
      self.prefix_re_matcher('MisformedToken'),
      self.prefix_re_matcher('MixtureModelError'),
      self.prefix_re_matcher('MotifsParseError'),
      self.prefix_re_matcher('MotifsSpaceError'),
      self.prefix_re_matcher('NoFastqDataException'),
      self.prefix_re_matcher('NotFullyVisible'),
      self.prefix_re_matcher('NullAxisMatrixError'),
      self.prefix_re_matcher('OsVersionError'),
      self.prefix_re_matcher('PreflightException'),
      self.prefix_re_matcher('ReadError'),
      self.prefix_re_matcher('ReanalyzeError'),
      self.prefix_re_matcher('RscriptError'),
      self.prefix_re_matcher('SimpleGoodTuringError'),
      self.prefix_re_matcher('SingleIndexFlowcellException'),
      self.prefix_re_matcher('UserInputError'),
      self.prefix_re_matcher('VDJReferenceConstructionError'),
      self.prefix_re_matcher('WhitelistNotFound'),
      // NumPy exceptions
      { exp: '^numpy\\._core._exceptions\\._ArrayMemoryError:', key: 'NumPy _ArrayMemoryError' },
      { exp: '^numpy\\._core._exceptions\\.', key: 'NumPy Error' },
      // Other Python exceptions
      { exp: '^\\S+Error:', key: 'Python Error' },
      { exp: '^\\S+Exception:', key: 'Python Exception' },
      // Rust assertion failed
      { exp: 'assertion `.*` failed', key: 'Rust assertion failed' },
      // Rust panic in Martian adapter
      { exp: '^stage failed unexpectedly', key: 'Rust panicked' },
      // Go signal error messages, either `Caught signal terminated` or `signal: terminated`
      { exp: 'signal:? hangup', key: 'SIGHUP' },  // 1
      { exp: 'signal:? interrupt', key: 'SIGINT' },  // 2
      { exp: 'signal:? quit', key: 'SIGQUIT' },  // 3
      { exp: 'signal:? illegal instruction', key: 'SIGILL' },  // 4
      { exp: 'signal:? trace/breakpoint trap', key: 'SIGTRAP' },  // 5
      { exp: 'signal:? aborted', key: 'SIGABRT' },  // 6
      { exp: 'signal:? bus error', key: 'SIGBUS' },  // 7
      { exp: 'signal:? floating point exception', key: 'SIGFPE' },  // 8
      { exp: 'signal:? killed', key: 'SIGKILL' },  // 9
      { exp: 'signal:? user defined signal 1', key: 'SIGUSR1' },  // 10
      { exp: 'signal:? segmentation fault', key: 'SIGSEGV' },  // 11
      { exp: 'signal:? user defined signal 2', key: 'SIGUSR2' },  // 12
      { exp: 'signal:? broken pipe', key: 'SIGPIPE' },  // 13
      { exp: 'signal:? alarm clock', key: 'SIGALRM' },  // 14
      { exp: 'signal:? terminated', key: 'SIGTERM' },  // 15
      { exp: 'signal:? stack fault', key: 'SIGSTKFLT' },  // 16
      { exp: 'signal:? child exited', key: 'SIGCHLD' },  // 17
      { exp: 'signal:? continued', key: 'SIGCONT' },  // 18
      { exp: 'signal:? stopped \\(signal\\)', key: 'SIGSTOP' },  // 19
      { exp: 'signal:? stopped', key: 'SIGTSTP' },  // 20
      { exp: 'signal:? stopped \\(tty input\\)', key: 'SIGTTIN' },  // 21
      { exp: 'signal:? stopped \\(tty output\\)', key: 'SIGTTOU' },  // 22
      { exp: 'signal:? urgent I/O condition', key: 'SIGURG' },  // 23
      { exp: 'signal:? CPU time limit exceeded', key: 'SIGXCPU' },  // 24
      { exp: 'signal:? file size limit exceeded', key: 'SIGXFSZ' },  // 25
      { exp: 'signal:? virtual timer expired', key: 'SIGVTALRM' },  // 26
      { exp: 'signal:? profiling timer expired', key: 'SIGPROF' },  // 27
      { exp: 'signal:? window changed', key: 'SIGWINCH' },  // 28
      { exp: 'signal:? I/O possible', key: 'SIGPOLL' },  // 29
      { exp: 'signal:? power failure', key: 'SIGPWR' },  // 30
      { exp: 'signal:? bad system call', key: 'SIGSYS' },  // 31
      // Other errors
      { exp: '^Error:', key: 'Error' },
    ],

  // Presets shared by all products.
  presets: {
    filesystem: {
      match: [
        {
          exp: '\\b(?:(?:9p)|(?:lxfs))\\b',
          key: 'wsl',
        },
        {
          exp: '\\b(?:(?:ext[34])|(?:(?:btr|x|z)fs))\\b',
          key: 'local',
        },
        {
          exp: '\\bnfs',
        },
        {
          exp: 'lustre',
        },
        {
          exp: 'overlay',
        },
        {
          exp: 'squashfs',
        },
        {
          exp: 'panfs',
        },
        {
          exp: 'gpfs',
        },
        {
          exp: 'beegfs',
        },
        {
          exp: 'cifs',
        },
        {
          exp: 'ceph',
        },
        {
          exp: 'fat\\b',
        },
        {
          exp: '\\bfuse',
        },
      ],
    },
    pca_v1_beta_coarse_cell_types: {
      match: [
        {
          exp: '(?i)^Low UMI Barcode$',
          key: 'low-umi-barcode',
        },
        {
          exp: '(?i)^abnormal cell$',
          key: 'abnormal-cell',
        },
        {
          exp: '(?i)^adipocyte$',
          key: 'adipocyte',
        },
        {
          exp: '(?i)^B cell$',
          key: 'b-cell',
        },
        {
          exp: '(?i)^cardiac muscle cell$',
          key: 'cardiac-muscle-cell',
        },
        {
          exp: '(?i)^cardiocyte$',
          key: 'cardiocyte',
        },
        {
          exp: '(?i)^ciliated cell$',
          key: 'ciliated-cell',
        },
        {
          exp: '(?i)^columnar/cuboidal epithelial cell$',
          key: 'columnar-cuboidal-epithelial-cell',
        },
        {
          exp: '(?i)^connective tissue cell$',
          key: 'connective-tissue-cell',
        },
        {
          exp: '(?i)^contractile cell$',
          key: 'contractile-cell',
        },
        {
          exp: '(?i)^endothelial cell$',
          key: 'endothelial-cell',
        },
        {
          exp: '(?i)^extraembryonic cell$',
          key: 'extraembryonic-cell',
        },
        {
          exp: '(?i)^fibroblast$',
          key: 'fibroblast',
        },
        {
          exp: '(?i)^GABAergic neuron$',
          key: 'gabaergic-neuron',
        },
        {
          exp: '(?i)^glandular epithelial cell$',
          key: 'glandular-epithelial-cell',
        },
        {
          exp: '(?i)^glial cell$',
          key: 'glial-cell',
        },
        {
          exp: '(?i)^glutamatergic neuron$',
          key: 'glutamatergic-neuron',
        },
        {
          exp: '(?i)^granulocyte$',
          key: 'granulocyte',
        },
        {
          exp: '(?i)^hematopoietic cell$',
          key: 'hematopoietic-cell',
        },
        {
          exp: '(?i)^hematopoietic precursor cell$',
          key: 'hematopoietic-precursor-cell',
        },
        {
          exp: '(?i)^lymphocyte$',
          key: 'lymphocyte',
        },
        {
          exp: '(?i)^macrophage$',
          key: 'macrophage',
        },
        {
          exp: '(?i)^mast cell$',
          key: 'mast-cell',
        },
        {
          exp: '(?i)^mesenchymal cell$',
          key: 'mesenchymal-cell',
        },
        {
          exp: '(?i)^monocyte$',
          key: 'monocyte',
        },
        {
          exp: '(?i)^mononuclear cell$',
          key: 'mononuclear-cell',
        },
        {
          exp: '(?i)^mural cell$',
          key: 'mural-cell',
        },
        {
          exp: '(?i)^myeloid leukocyte$',
          key: 'myeloid-leukocyte',
        },
        {
          exp: '(?i)^leukocyte$',
          key: 'leukocyte',
        },
        {
          exp: '(?i)^neural cell$',
          key: 'neural-cell',
        },
        {
          exp: '(?i)^neuron associated cell$',
          key: 'neuron-associated-cell',
        },
        {
          exp: '(?i)^neuron$',
          key: 'neuron',
        },
        {
          exp: '(?i)^perivascular cell$',
          key: 'perivascular-cell',
        },
        {
          exp: '(?i)^precursor cell$',
          key: 'precursor-cell',
        },
        {
          exp: '(?i)^primordial germ cell$',
          key: 'primordial-germ-cell',
        },
        {
          exp: '(?i)^germ cell$',
          key: 'germ-cell',
        },
        {
          exp: '(?i)^salivary gland cell$',
          key: 'salivary-gland-cell',
        },
        {
          exp: '(?i)^sebaceous gland cell$',
          key: 'sebaceous-gland-cell',
        },
        {
          exp: '(?i)^secretory cell$',
          key: 'secretory-cell',
        },
        {
          exp: '(?i)^serous secreting cell$',
          key: 'serous-secreting-cell',
        },
        {
          exp: '(?i)^smooth muscle cell$',
          key: 'smooth-muscle-cell',
        },
        {
          exp: '(?i)^muscle cell$',
          key: 'muscle-cell',
        },
        {
          exp: '(?i)^squamous epithelial cell$',
          key: 'squamous-epithelial-cell',
        },
        {
          exp: '(?i)^epithelial cell$',
          key: 'epithelial-cell',
        },
        {
          exp: '(?i)^stem cell$',
          key: 'stem-cell',
        },
        {
          exp: '(?i)^T cell$',
          key: 't-cell',
        },
        {
          exp: '(?i)^transit amplifying cell$',
          key: 'transit-amplifying-cell',
        },
        {
          exp: '(?i)^cell$',
          key: 'cell',
        },
      ],
    },
    pca_v1_beta_fine_cell_types: {
      match: [
        {
          exp: '(?i)^Low UMI Barcode$',
          key: 'low-umi-barcode',
        },
        {
          exp: '(?i)^A2 amacrine cell$',
          key: 'a2-amacrine-cell',
        },
        {
          exp: '(?i)^B-1 B cell$',
          key: 'b-1-b-cell',
        },
        {
          exp: '(?i)^B-1a B cell$',
          key: 'b-1a-b-cell',
        },
        {
          exp: '(?i)^B-1b B cell$',
          key: 'b-1b-b-cell',
        },
        {
          exp: '(?i)^B-2 B cell$',
          key: 'b-2-b-cell',
        },
        {
          exp: '(?i)^BEST4\\+ intestinal epithelial cell, human$',
          key: 'best4-intestinal-epithelial-cell-human',
        },
        {
          exp: '(?i)^Bergmann glial cell$',
          key: 'bergmann-glial-cell',
        },
        {
          exp: '(?i)^CD14-low, CD16-positive monocyte$',
          key: 'cd14-low-cd16-positive-monocyte',
        },
        {
          exp: '(?i)^CD14-positive monocyte$',
          key: 'cd14-positive-monocyte',
        },
        {
          exp: '(?i)^CD14-positive, CD16-negative classical monocyte$',
          key: 'cd14-positive-cd16-negative-classical-monocyte',
        },
        {
          exp: '(?i)^CD14-positive, CD16-positive monocyte$',
          key: 'cd14-positive-cd16-positive-monocyte',
        },
        {
          exp: '(?i)^CD141-positive myeloid dendritic cell$',
          key: 'cd141-positive-myeloid-dendritic-cell',
        },
        {
          exp: '(?i)^CD16-negative, CD56-bright natural killer cell, human$',
          key: 'cd16-negative-cd56-bright-natural-killer-cell-human',
        },
        {
          exp: '(?i)^CD16-positive, CD56-dim natural killer cell, human$',
          key: 'cd16-positive-cd56-dim-natural-killer-cell-human',
        },
        {
          exp: '(?i)^CD1c-positive myeloid dendritic cell$',
          key: 'cd1c-positive-myeloid-dendritic-cell',
        },
        {
          exp: '(?i)^CD34-positive, CD38-negative hematopoietic stem cell$',
          key: 'cd34-positive-cd38-negative-hematopoietic-stem-cell',
        },
        {
          exp: '(?i)^CD34-positive, CD56-positive, CD117-positive common innate lymphoid precursor, human$',
          key: 'cd34-positive-cd56-positive-cd117-positive-common-innate-lymphoid-precursor-human',
        },
        {
          exp: '(?i)^CD4-positive helper T cell$',
          key: 'cd4-positive-helper-t-cell',
        },
        {
          exp: '(?i)^CD4-positive, CD25-positive, alpha-beta regulatory T cell$',
          key: 'cd4-positive-cd25-positive-alpha-beta-regulatory-t-cell',
        },
        {
          exp: '(?i)^CD4-positive, alpha-beta cytotoxic T cell$',
          key: 'cd4-positive-alpha-beta-cytotoxic-t-cell',
        },
        {
          exp: '(?i)^CD4-positive, alpha-beta memory T cell$',
          key: 'cd4-positive-alpha-beta-memory-t-cell',
        },
        {
          exp: '(?i)^CD4-positive, alpha-beta thymocyte$',
          key: 'cd4-positive-alpha-beta-thymocyte',
        },
        {
          exp: '(?i)^CD56-positive, CD161-positive immature natural killer cell, human$',
          key: 'cd56-positive-cd161-positive-immature-natural-killer-cell-human',
        },
        {
          exp: '(?i)^CD8-alpha alpha positive, gamma-delta intraepithelial T cell$',
          key: 'cd8-alpha-alpha-positive-gamma-delta-intraepithelial-t-cell',
        },
        {
          exp: '(?i)^CD8-alpha-alpha-positive, alpha-beta intraepithelial T cell$',
          key: 'cd8-alpha-alpha-positive-alpha-beta-intraepithelial-t-cell',
        },
        {
          exp: '(?i)^CD8-alpha-beta-positive, alpha-beta intraepithelial T cell$',
          key: 'cd8-alpha-beta-positive-alpha-beta-intraepithelial-t-cell',
        },
        {
          exp: '(?i)^CD8-positive, alpha-beta cytokine secreting effector T cell$',
          key: 'cd8-positive-alpha-beta-cytokine-secreting-effector-t-cell',
        },
        {
          exp: '(?i)^CD8-positive, alpha-beta cytotoxic T cell$',
          key: 'cd8-positive-alpha-beta-cytotoxic-t-cell',
        },
        {
          exp: '(?i)^CD8-positive, alpha-beta memory T cell, CD45RO-positive$',
          key: 'cd8-positive-alpha-beta-memory-t-cell-cd45ro-positive',
        },
        {
          exp: '(?i)^CD8-positive, alpha-beta memory T cell$',
          key: 'cd8-positive-alpha-beta-memory-t-cell',
        },
        {
          exp: '(?i)^CD8-positive, alpha-beta thymocyte$',
          key: 'cd8-positive-alpha-beta-thymocyte',
        },
        {
          exp: '(?i)^Cajal-Retzius cell$',
          key: 'cajal-retzius-cell',
        },
        {
          exp: '(?i)^DN1 thymic pro-T cell$',
          key: 'dn1-thymic-pro-t-cell',
        },
        {
          exp: '(?i)^DN3 thymocyte$',
          key: 'dn3-thymocyte',
        },
        {
          exp: '(?i)^DN4 thymocyte$',
          key: 'dn4-thymocyte',
        },
        {
          exp: '(?i)^GABAergic amacrine cell$',
          key: 'gabaergic-amacrine-cell',
        },
        {
          exp: '(?i)^GABAergic neuron$',
          key: 'gabaergic-neuron',
        },
        {
          exp: '(?i)^GIP cell$',
          key: 'gip-cell',
        },
        {
          exp: '(?i)^H1 horizontal cell$',
          key: 'h1-horizontal-cell',
        },
        {
          exp: '(?i)^H2 horizontal cell$',
          key: 'h2-horizontal-cell',
        },
        {
          exp: '(?i)^Hofbauer cell$',
          key: 'hofbauer-cell',
        },
        {
          exp: '(?i)^ILC1, human$',
          key: 'ilc1-human',
        },
        {
          exp: '(?i)^IgA plasma cell$',
          key: 'iga-plasma-cell',
        },
        {
          exp: '(?i)^IgA plasmablast$',
          key: 'iga-plasmablast',
        },
        {
          exp: '(?i)^IgG memory B cell$',
          key: 'igg-memory-b-cell',
        },
        {
          exp: '(?i)^IgG plasma cell$',
          key: 'igg-plasma-cell',
        },
        {
          exp: '(?i)^IgG plasmablast$',
          key: 'igg-plasmablast',
        },
        {
          exp: '(?i)^IgG-negative class switched memory B cell$',
          key: 'igg-negative-class-switched-memory-b-cell',
        },
        {
          exp: '(?i)^IgM plasma cell$',
          key: 'igm-plasma-cell',
        },
        {
          exp: '(?i)^Kupffer cell$',
          key: 'kupffer-cell',
        },
        {
          exp: '(?i)^L2/3 intratelencephalic projecting glutamatergic neuron of the primary motor cortex$',
          key: 'l2-3-intratelencephalic-projecting-glutamatergic-neuron-of-the-primary-motor-cortex',
        },
        {
          exp: '(?i)^L2/3 intratelencephalic projecting glutamatergic neuron$',
          key: 'l2-3-intratelencephalic-projecting-glutamatergic-neuron',
        },
        {
          exp: '(?i)^L2/3-6 intratelencephalic projecting glutamatergic neuron$',
          key: 'l2-3-6-intratelencephalic-projecting-glutamatergic-neuron',
        },
        {
          exp: '(?i)^L4 intratelencephalic projecting glutamatergic neuron$',
          key: 'l4-intratelencephalic-projecting-glutamatergic-neuron',
        },
        {
          exp: '(?i)^L4/5 intratelencephalic projecting glutamatergic neuron of the primary motor cortex$',
          key: 'l4-5-intratelencephalic-projecting-glutamatergic-neuron-of-the-primary-motor-cortex',
        },
        {
          exp: '(?i)^L5 extratelencephalic projecting glutamatergic cortical neuron$',
          key: 'l5-extratelencephalic-projecting-glutamatergic-cortical-neuron',
        },
        {
          exp: '(?i)^L5/6 near-projecting glutamatergic neuron of the primary motor cortex$',
          key: 'l5-6-near-projecting-glutamatergic-neuron-of-the-primary-motor-cortex',
        },
        {
          exp: '(?i)^L5/6 near-projecting glutamatergic neuron$',
          key: 'l5-6-near-projecting-glutamatergic-neuron',
        },
        {
          exp: '(?i)^L6 corticothalamic-projecting glutamatergic cortical neuron$',
          key: 'l6-corticothalamic-projecting-glutamatergic-cortical-neuron',
        },
        {
          exp: '(?i)^L6 intratelencephalic projecting glutamatergic neuron of the primary motor cortex$',
          key: 'l6-intratelencephalic-projecting-glutamatergic-neuron-of-the-primary-motor-cortex',
        },
        {
          exp: '(?i)^L6b glutamatergic cortical neuron$',
          key: 'l6b-glutamatergic-cortical-neuron',
        },
        {
          exp: '(?i)^Leydig cell$',
          key: 'leydig-cell',
        },
        {
          exp: '(?i)^M cell of gut$',
          key: 'm-cell-of-gut',
        },
        {
          exp: '(?i)^MHC-II-positive classical monocyte$',
          key: 'mhc-ii-positive-classical-monocyte',
        },
        {
          exp: '(?i)^Merkel cell$',
          key: 'merkel-cell',
        },
        {
          exp: '(?i)^Mueller cell$',
          key: 'mueller-cell',
        },
        {
          exp: '(?i)^NKp44-negative group 3 innate lymphoid cell, human$',
          key: 'nkp44-negative-group-3-innate-lymphoid-cell-human',
        },
        {
          exp: '(?i)^NKp44-positive group 3 innate lymphoid cell, human$',
          key: 'nkp44-positive-group-3-innate-lymphoid-cell-human',
        },
        {
          exp: '(?i)^NKp46-positive innate lymphoid cell, human$',
          key: 'nkp46-positive-innate-lymphoid-cell-human',
        },
        {
          exp: '(?i)^OFF retinal ganglion cell$',
          key: 'off-retinal-ganglion-cell',
        },
        {
          exp: '(?i)^OFF-bipolar cell$',
          key: 'off-bipolar-cell',
        },
        {
          exp: '(?i)^OFFx cell$',
          key: 'offx-cell',
        },
        {
          exp: '(?i)^ON retinal ganglion cell$',
          key: 'on-retinal-ganglion-cell',
        },
        {
          exp: '(?i)^ON-bipolar cell$',
          key: 'on-bipolar-cell',
        },
        {
          exp: '(?i)^ON-blue cone bipolar cell$',
          key: 'on-blue-cone-bipolar-cell',
        },
        {
          exp: '(?i)^P/D1 enteroendocrine cell$',
          key: 'p-d1-enteroendocrine-cell',
        },
        {
          exp: '(?i)^Purkinje cell$',
          key: 'purkinje-cell',
        },
        {
          exp: '(?i)^S cone cell$',
          key: 's-cone-cell',
        },
        {
          exp: '(?i)^Schwann cell precursor$',
          key: 'schwann-cell-precursor',
        },
        {
          exp: '(?i)^Sertoli cell$',
          key: 'sertoli-cell',
        },
        {
          exp: '(?i)^T follicular helper cell$',
          key: 't-follicular-helper-cell',
        },
        {
          exp: '(?i)^T follicular regulatory cell$',
          key: 't-follicular-regulatory-cell',
        },
        {
          exp: '(?i)^T-helper 1 cell$',
          key: 't-helper-1-cell',
        },
        {
          exp: '(?i)^T-helper 17 cell$',
          key: 't-helper-17-cell',
        },
        {
          exp: '(?i)^T-helper 2 cell$',
          key: 't-helper-2-cell',
        },
        {
          exp: '(?i)^T-helper 22 cell$',
          key: 't-helper-22-cell',
        },
        {
          exp: '(?i)^Tc1 cell$',
          key: 'tc1-cell',
        },
        {
          exp: '(?i)^abnormal cell$',
          key: 'abnormal-cell',
        },
        {
          exp: '(?i)^acinar cell of salivary gland$',
          key: 'acinar-cell-of-salivary-gland',
        },
        {
          exp: '(?i)^activated CD4-positive, alpha-beta T cell, human$',
          key: 'activated-cd4-positive-alpha-beta-t-cell-human',
        },
        {
          exp: '(?i)^activated CD4-positive, alpha-beta T cell$',
          key: 'activated-cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^activated CD8-positive, alpha-beta T cell, human$',
          key: 'activated-cd8-positive-alpha-beta-t-cell-human',
        },
        {
          exp: '(?i)^activated CD8-positive, alpha-beta T cell$',
          key: 'activated-cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^activated type II NK T cell$',
          key: 'activated-type-ii-nk-t-cell',
        },
        {
          exp: '(?i)^adipocyte of breast$',
          key: 'adipocyte-of-breast',
        },
        {
          exp: '(?i)^adipocyte of epicardial fat of left ventricle$',
          key: 'adipocyte-of-epicardial-fat-of-left-ventricle',
        },
        {
          exp: '(?i)^adventitial cell$',
          key: 'adventitial-cell',
        },
        {
          exp: '(?i)^airway submucosal gland duct basal cell$',
          key: 'airway-submucosal-gland-duct-basal-cell',
        },
        {
          exp: '(?i)^alveolar type 1 fibroblast cell$',
          key: 'alveolar-type-1-fibroblast-cell',
        },
        {
          exp: '(?i)^alveolar type 2 fibroblast cell$',
          key: 'alveolar-type-2-fibroblast-cell',
        },
        {
          exp: '(?i)^anterior lens cell$',
          key: 'anterior-lens-cell',
        },
        {
          exp: '(?i)^aortic endothelial cell$',
          key: 'aortic-endothelial-cell',
        },
        {
          exp: '(?i)^aortic smooth muscle cell$',
          key: 'aortic-smooth-muscle-cell',
        },
        {
          exp: '(?i)^astrocyte of the cerebellum$',
          key: 'astrocyte-of-the-cerebellum',
        },
        {
          exp: '(?i)^astrocyte of the cerebral cortex$',
          key: 'astrocyte-of-the-cerebral-cortex',
        },
        {
          exp: '(?i)^basal cell of epidermis$',
          key: 'basal-cell-of-epidermis',
        },
        {
          exp: '(?i)^basal cell of epithelium of trachea$',
          key: 'basal-cell-of-epithelium-of-trachea',
        },
        {
          exp: '(?i)^basal cell of prostate epithelium$',
          key: 'basal-cell-of-prostate-epithelium',
        },
        {
          exp: '(?i)^basal epithelial cell of prostatic duct$',
          key: 'basal-epithelial-cell-of-prostatic-duct',
        },
        {
          exp: '(?i)^basal epithelial cell of tracheobronchial tree$',
          key: 'basal-epithelial-cell-of-tracheobronchial-tree',
        },
        {
          exp: '(?i)^basophil mast progenitor cell$',
          key: 'basophil-mast-progenitor-cell',
        },
        {
          exp: '(?i)^basophil$',
          key: 'basophil',
        },
        {
          exp: '(?i)^bladder cell$',
          key: 'bladder-cell',
        },
        {
          exp: '(?i)^bladder urothelial cell$',
          key: 'bladder-urothelial-cell',
        },
        {
          exp: '(?i)^blood vessel smooth muscle cell$',
          key: 'blood-vessel-smooth-muscle-cell',
        },
        {
          exp: '(?i)^bone marrow cell$',
          key: 'bone-marrow-cell',
        },
        {
          exp: '(?i)^brain pericyte$',
          key: 'brain-pericyte',
        },
        {
          exp: '(?i)^brain vascular cell$',
          key: 'brain-vascular-cell',
        },
        {
          exp: '(?i)^bronchial epithelial cell$',
          key: 'bronchial-epithelial-cell',
        },
        {
          exp: '(?i)^bronchus fibroblast of lung$',
          key: 'bronchus-fibroblast-of-lung',
        },
        {
          exp: '(?i)^brush cell of epithelium proper of large intestine$',
          key: 'brush-cell-of-epithelium-proper-of-large-intestine',
        },
        {
          exp: '(?i)^brush cell of trachea$',
          key: 'brush-cell-of-trachea',
        },
        {
          exp: '(?i)^brush cell of trachebronchial tree$',
          key: 'brush-cell-of-trachebronchial-tree',
        },
        {
          exp: '(?i)^cardiac endothelial cell$',
          key: 'cardiac-endothelial-cell',
        },
        {
          exp: '(?i)^cardiac mesenchymal cell$',
          key: 'cardiac-mesenchymal-cell',
        },
        {
          exp: '(?i)^cardiac muscle myoblast$',
          key: 'cardiac-muscle-myoblast',
        },
        {
          exp: '(?i)^cardiac neuron$',
          key: 'cardiac-neuron',
        },
        {
          exp: '(?i)^cardiac valve cell$',
          key: 'cardiac-valve-cell',
        },
        {
          exp: '(?i)^cardiocyte$',
          key: 'cardiocyte',
        },
        {
          exp: '(?i)^caudal ganglionic eminence derived GABAergic cortical interneuron$',
          key: 'caudal-ganglionic-eminence-derived-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^cell of skeletal muscle$',
          key: 'cell-of-skeletal-muscle',
        },
        {
          exp: '(?i)^central memory CD4-positive, alpha-beta T cell$',
          key: 'central-memory-cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^central memory CD8-positive, alpha-beta T cell$',
          key: 'central-memory-cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^central nervous system macrophage$',
          key: 'central-nervous-system-macrophage',
        },
        {
          exp: '(?i)^centrilobular region hepatocyte$',
          key: 'centrilobular-region-hepatocyte',
        },
        {
          exp: '(?i)^centroblast$',
          key: 'centroblast',
        },
        {
          exp: '(?i)^centrocyte$',
          key: 'centrocyte',
        },
        {
          exp: '(?i)^cerebellar Golgi cell$',
          key: 'cerebellar-golgi-cell',
        },
        {
          exp: '(?i)^cerebellar granule cell precursor$',
          key: 'cerebellar-granule-cell-precursor',
        },
        {
          exp: '(?i)^cerebellar granule cell$',
          key: 'cerebellar-granule-cell',
        },
        {
          exp: '(?i)^cerebral cortex GABAergic interneuron$',
          key: 'cerebral-cortex-gabaergic-interneuron',
        },
        {
          exp: '(?i)^cerebral cortex endothelial cell$',
          key: 'cerebral-cortex-endothelial-cell',
        },
        {
          exp: '(?i)^cerebral cortex neuron$',
          key: 'cerebral-cortex-neuron',
        },
        {
          exp: '(?i)^chandelier cell$',
          key: 'chandelier-cell',
        },
        {
          exp: '(?i)^chandelier pvalb GABAergic cortical interneuron$',
          key: 'chandelier-pvalb-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^chondrocyte$',
          key: 'chondrocyte',
        },
        {
          exp: '(?i)^chorionic trophoblast cell$',
          key: 'chorionic-trophoblast-cell',
        },
        {
          exp: '(?i)^choroid plexus epithelial cell$',
          key: 'choroid-plexus-epithelial-cell',
        },
        {
          exp: '(?i)^chromaffin cell$',
          key: 'chromaffin-cell',
        },
        {
          exp: '(?i)^ciliary muscle cell$',
          key: 'ciliary-muscle-cell',
        },
        {
          exp: '(?i)^ciliated cell of the bronchus$',
          key: 'ciliated-cell-of-the-bronchus',
        },
        {
          exp: '(?i)^ciliated columnar cell of tracheobronchial tree$',
          key: 'ciliated-columnar-cell-of-tracheobronchial-tree',
        },
        {
          exp: '(?i)^class switched memory B cell$',
          key: 'class-switched-memory-b-cell',
        },
        {
          exp: '(?i)^club cell$',
          key: 'club-cell',
        },
        {
          exp: '(?i)^collagen secreting cell$',
          key: 'collagen-secreting-cell',
        },
        {
          exp: '(?i)^colon epithelial cell$',
          key: 'colon-epithelial-cell',
        },
        {
          exp: '(?i)^colon goblet cell$',
          key: 'colon-goblet-cell',
        },
        {
          exp: '(?i)^colon macrophage$',
          key: 'colon-macrophage',
        },
        {
          exp: '(?i)^columnar/cuboidal epithelial cell$',
          key: 'columnar-cuboidal-epithelial-cell',
        },
        {
          exp: '(?i)^common dendritic progenitor$',
          key: 'common-dendritic-progenitor',
        },
        {
          exp: '(?i)^common lymphoid progenitor$',
          key: 'common-lymphoid-progenitor',
        },
        {
          exp: '(?i)^common myeloid progenitor$',
          key: 'common-myeloid-progenitor',
        },
        {
          exp: '(?i)^cone retinal bipolar cell$',
          key: 'cone-retinal-bipolar-cell',
        },
        {
          exp: '(?i)^conjunctival epithelial cell$',
          key: 'conjunctival-epithelial-cell',
        },
        {
          exp: '(?i)^connective tissue cell$',
          key: 'connective-tissue-cell',
        },
        {
          exp: '(?i)^contractile cell$',
          key: 'contractile-cell',
        },
        {
          exp: '(?i)^cord blood hematopoietic stem cell$',
          key: 'cord-blood-hematopoietic-stem-cell',
        },
        {
          exp: '(?i)^corneal endothelial cell$',
          key: 'corneal-endothelial-cell',
        },
        {
          exp: '(?i)^corneal epithelial cell$',
          key: 'corneal-epithelial-cell',
        },
        {
          exp: '(?i)^cortical cell of adrenal gland$',
          key: 'cortical-cell-of-adrenal-gland',
        },
        {
          exp: '(?i)^cortical thymic epithelial cell$',
          key: 'cortical-thymic-epithelial-cell',
        },
        {
          exp: '(?i)^corticothalamic-projecting glutamatergic cortical neuron$',
          key: 'corticothalamic-projecting-glutamatergic-cortical-neuron',
        },
        {
          exp: '(?i)^cytotoxic T cell$',
          key: 'cytotoxic-t-cell',
        },
        {
          exp: '(?i)^decidual cell$',
          key: 'decidual-cell',
        },
        {
          exp: '(?i)^decidual natural killer cell, human$',
          key: 'decidual-natural-killer-cell-human',
        },
        {
          exp: '(?i)^dermis microvascular lymphatic vessel endothelial cell$',
          key: 'dermis-microvascular-lymphatic-vessel-endothelial-cell',
        },
        {
          exp: '(?i)^deuterosomal cell$',
          key: 'deuterosomal-cell',
        },
        {
          exp: '(?i)^differentiation-committed oligodendrocyte precursor$',
          key: 'differentiation-committed-oligodendrocyte-precursor',
        },
        {
          exp: '(?i)^diffuse bipolar 1 cell$',
          key: 'diffuse-bipolar-1-cell',
        },
        {
          exp: '(?i)^diffuse bipolar 2 cell$',
          key: 'diffuse-bipolar-2-cell',
        },
        {
          exp: '(?i)^diffuse bipolar 3a cell$',
          key: 'diffuse-bipolar-3a-cell',
        },
        {
          exp: '(?i)^diffuse bipolar 3b cell$',
          key: 'diffuse-bipolar-3b-cell',
        },
        {
          exp: '(?i)^diffuse bipolar 4 cell$',
          key: 'diffuse-bipolar-4-cell',
        },
        {
          exp: '(?i)^diffuse bipolar 6 cell$',
          key: 'diffuse-bipolar-6-cell',
        },
        {
          exp: '(?i)^dopaminergic neuron$',
          key: 'dopaminergic-neuron',
        },
        {
          exp: '(?i)^double negative T regulatory cell$',
          key: 'double-negative-t-regulatory-cell',
        },
        {
          exp: '(?i)^double negative thymocyte$',
          key: 'double-negative-thymocyte',
        },
        {
          exp: '(?i)^double-positive, alpha-beta thymocyte$',
          key: 'double-positive-alpha-beta-thymocyte',
        },
        {
          exp: '(?i)^duodenum glandular cell$',
          key: 'duodenum-glandular-cell',
        },
        {
          exp: '(?i)^early T lineage precursor$',
          key: 'early-t-lineage-precursor',
        },
        {
          exp: '(?i)^early lymphoid progenitor$',
          key: 'early-lymphoid-progenitor',
        },
        {
          exp: '(?i)^early pro-B cell$',
          key: 'early-pro-b-cell',
        },
        {
          exp: '(?i)^early promyelocyte$',
          key: 'early-promyelocyte',
        },
        {
          exp: '(?i)^effector CD4-positive, alpha-beta T cell$',
          key: 'effector-cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^effector CD8-positive, alpha-beta T cell$',
          key: 'effector-cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^effector memory CD4-positive, alpha-beta T cell$',
          key: 'effector-memory-cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^effector memory CD8-positive, alpha-beta T cell, terminally differentiated$',
          key: 'effector-memory-cd8-positive-alpha-beta-t-cell-terminally-differentiated',
        },
        {
          exp: '(?i)^effector memory CD8-positive, alpha-beta T cell$',
          key: 'effector-memory-cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^elicited macrophage$',
          key: 'elicited-macrophage',
        },
        {
          exp: '(?i)^embryonic fibroblast$',
          key: 'embryonic-fibroblast',
        },
        {
          exp: '(?i)^embryonic stem cell$',
          key: 'embryonic-stem-cell',
        },
        {
          exp: '(?i)^endocardial cell$',
          key: 'endocardial-cell',
        },
        {
          exp: '(?i)^endosteal cell$',
          key: 'endosteal-cell',
        },
        {
          exp: '(?i)^endothelial cell of artery$',
          key: 'endothelial-cell-of-artery',
        },
        {
          exp: '(?i)^endothelial cell of coronary artery$',
          key: 'endothelial-cell-of-coronary-artery',
        },
        {
          exp: '(?i)^endothelial cell of hepatic sinusoid$',
          key: 'endothelial-cell-of-hepatic-sinusoid',
        },
        {
          exp: '(?i)^endothelial cell of lymphatic vessel$',
          key: 'endothelial-cell-of-lymphatic-vessel',
        },
        {
          exp: '(?i)^endothelial cell of pericentral hepatic sinusoid$',
          key: 'endothelial-cell-of-pericentral-hepatic-sinusoid',
        },
        {
          exp: '(?i)^endothelial cell of periportal hepatic sinusoid$',
          key: 'endothelial-cell-of-periportal-hepatic-sinusoid',
        },
        {
          exp: '(?i)^endothelial cell of placenta$',
          key: 'endothelial-cell-of-placenta',
        },
        {
          exp: '(?i)^endothelial cell of sinusoid$',
          key: 'endothelial-cell-of-sinusoid',
        },
        {
          exp: '(?i)^endothelial cell of uterus$',
          key: 'endothelial-cell-of-uterus',
        },
        {
          exp: '(?i)^endothelial cell of vascular tree$',
          key: 'endothelial-cell-of-vascular-tree',
        },
        {
          exp: '(?i)^endothelial cell of venule$',
          key: 'endothelial-cell-of-venule',
        },
        {
          exp: '(?i)^endothelial stalk cell$',
          key: 'endothelial-stalk-cell',
        },
        {
          exp: '(?i)^endothelial tip cell$',
          key: 'endothelial-tip-cell',
        },
        {
          exp: '(?i)^enteric neuron$',
          key: 'enteric-neuron',
        },
        {
          exp: '(?i)^enteric smooth muscle cell$',
          key: 'enteric-smooth-muscle-cell',
        },
        {
          exp: '(?i)^enterocyte of colon$',
          key: 'enterocyte-of-colon',
        },
        {
          exp: '(?i)^enterocyte of epithelium of large intestine$',
          key: 'enterocyte-of-epithelium-of-large-intestine',
        },
        {
          exp: '(?i)^enterocyte of epithelium of small intestine$',
          key: 'enterocyte-of-epithelium-of-small-intestine',
        },
        {
          exp: '(?i)^enterocyte of epithelium proper of ileum$',
          key: 'enterocyte-of-epithelium-proper-of-ileum',
        },
        {
          exp: '(?i)^enterocyte$',
          key: 'enterocyte',
        },
        {
          exp: '(?i)^enteroendocrine cell of colon$',
          key: 'enteroendocrine-cell-of-colon',
        },
        {
          exp: '(?i)^enteroendocrine cell of small intestine$',
          key: 'enteroendocrine-cell-of-small-intestine',
        },
        {
          exp: '(?i)^enucleate erythrocyte$',
          key: 'enucleate-erythrocyte',
        },
        {
          exp: '(?i)^enucleated reticulocyte$',
          key: 'enucleated-reticulocyte',
        },
        {
          exp: '(?i)^eosinophil$',
          key: 'eosinophil',
        },
        {
          exp: '(?i)^ependymal cell$',
          key: 'ependymal-cell',
        },
        {
          exp: '(?i)^epicardial adipocyte$',
          key: 'epicardial-adipocyte',
        },
        {
          exp: '(?i)^epidermal Langerhans cell$',
          key: 'epidermal-langerhans-cell',
        },
        {
          exp: '(?i)^Langerhans cell$',
          key: 'langerhans-cell',
        },
        {
          exp: '(?i)^epidermal cell$',
          key: 'epidermal-cell',
        },
        {
          exp: '(?i)^epithelial cell of alveolus of lung$',
          key: 'epithelial-cell-of-alveolus-of-lung',
        },
        {
          exp: '(?i)^epithelial cell of esophagus$',
          key: 'epithelial-cell-of-esophagus',
        },
        {
          exp: '(?i)^epithelial cell of exocrine pancreas$',
          key: 'epithelial-cell-of-exocrine-pancreas',
        },
        {
          exp: '(?i)^epithelial cell of glomerular capsule$',
          key: 'epithelial-cell-of-glomerular-capsule',
        },
        {
          exp: '(?i)^epithelial cell of lacrimal sac$',
          key: 'epithelial-cell-of-lacrimal-sac',
        },
        {
          exp: '(?i)^epithelial cell of large intestine$',
          key: 'epithelial-cell-of-large-intestine',
        },
        {
          exp: '(?i)^epithelial cell of lower respiratory tract$',
          key: 'epithelial-cell-of-lower-respiratory-tract',
        },
        {
          exp: '(?i)^epithelial cell of lung$',
          key: 'epithelial-cell-of-lung',
        },
        {
          exp: '(?i)^epithelial cell of nephron$',
          key: 'epithelial-cell-of-nephron',
        },
        {
          exp: '(?i)^epithelial cell of prostate$',
          key: 'epithelial-cell-of-prostate',
        },
        {
          exp: '(?i)^epithelial cell of proximal tubule segment 1$',
          key: 'epithelial-cell-of-proximal-tubule-segment-1',
        },
        {
          exp: '(?i)^epithelial cell of proximal tubule segment 3$',
          key: 'epithelial-cell-of-proximal-tubule-segment-3',
        },
        {
          exp: '(?i)^epithelial cell of proximal tubule$',
          key: 'epithelial-cell-of-proximal-tubule',
        },
        {
          exp: '(?i)^epithelial cell of small intestine$',
          key: 'epithelial-cell-of-small-intestine',
        },
        {
          exp: '(?i)^epithelial cell of stratum germinativum of esophagus$',
          key: 'epithelial-cell-of-stratum-germinativum-of-esophagus',
        },
        {
          exp: '(?i)^epithelial cell of sweat gland$',
          key: 'epithelial-cell-of-sweat-gland',
        },
        {
          exp: '(?i)^epithelial cell of thymus$',
          key: 'epithelial-cell-of-thymus',
        },
        {
          exp: '(?i)^epithelial cell of urethra$',
          key: 'epithelial-cell-of-urethra',
        },
        {
          exp: '(?i)^epithelial cell of uterus$',
          key: 'epithelial-cell-of-uterus',
        },
        {
          exp: '(?i)^epithelial fate stem cell$',
          key: 'epithelial-fate-stem-cell',
        },
        {
          exp: '(?i)^erythrocyte$',
          key: 'erythrocyte',
        },
        {
          exp: '(?i)^erythroid lineage cell$',
          key: 'erythroid-lineage-cell',
        },
        {
          exp: '(?i)^erythroid progenitor cell, mammalian$',
          key: 'erythroid-progenitor-cell-mammalian',
        },
        {
          exp: '(?i)^eukaryotic cell$',
          key: 'eukaryotic-cell',
        },
        {
          exp: '(?i)^eurydendroid cell$',
          key: 'eurydendroid-cell',
        },
        {
          exp: '(?i)^exhausted T cell$',
          key: 'exhausted-t-cell',
        },
        {
          exp: '(?i)^extraembryonic cell$',
          key: 'extraembryonic-cell',
        },
        {
          exp: '(?i)^extravillous trophoblast$',
          key: 'extravillous-trophoblast',
        },
        {
          exp: '(?i)^eye photoreceptor cell$',
          key: 'eye-photoreceptor-cell',
        },
        {
          exp: '(?i)^fast muscle cell$',
          key: 'fast-muscle-cell',
        },
        {
          exp: '(?i)^female germ cell$',
          key: 'female-germ-cell',
        },
        {
          exp: '(?i)^fenestrated cell$',
          key: 'fenestrated-cell',
        },
        {
          exp: '(?i)^fetal cardiomyocyte$',
          key: 'fetal-cardiomyocyte',
        },
        {
          exp: '(?i)^fibro/adipogenic progenitor cell$',
          key: 'fibro-adipogenic-progenitor-cell',
        },
        {
          exp: '(?i)^fibroblast of breast$',
          key: 'fibroblast-of-breast',
        },
        {
          exp: '(?i)^fibroblast of cardiac tissue$',
          key: 'fibroblast-of-cardiac-tissue',
        },
        {
          exp: '(?i)^fibroblast of connective tissue of glandular part of prostate$',
          key: 'fibroblast-of-connective-tissue-of-glandular-part-of-prostate',
        },
        {
          exp: '(?i)^fibroblast of connective tissue of nonglandular part of prostate$',
          key: 'fibroblast-of-connective-tissue-of-nonglandular-part-of-prostate',
        },
        {
          exp: '(?i)^fibroblast of connective tissue of prostate$',
          key: 'fibroblast-of-connective-tissue-of-prostate',
        },
        {
          exp: '(?i)^fibroblast of lung$',
          key: 'fibroblast-of-lung',
        },
        {
          exp: '(?i)^fibroblast of lymphatic vessel$',
          key: 'fibroblast-of-lymphatic-vessel',
        },
        {
          exp: '(?i)^fibroblast of mammary gland$',
          key: 'fibroblast-of-mammary-gland',
        },
        {
          exp: '(?i)^fibrocyte$',
          key: 'fibrocyte',
        },
        {
          exp: '(?i)^flat midget bipolar cell$',
          key: 'flat-midget-bipolar-cell',
        },
        {
          exp: '(?i)^follicular B cell$',
          key: 'follicular-b-cell',
        },
        {
          exp: '(?i)^follicular dendritic cell$',
          key: 'follicular-dendritic-cell',
        },
        {
          exp: '(?i)^forebrain neuroblast$',
          key: 'forebrain-neuroblast',
        },
        {
          exp: '(?i)^forebrain radial glial cell$',
          key: 'forebrain-radial-glial-cell',
        },
        {
          exp: '(?i)^foveolar cell of stomach$',
          key: 'foveolar-cell-of-stomach',
        },
        {
          exp: '(?i)^fraction A pre-pro B cell$',
          key: 'fraction-a-pre-pro-b-cell',
        },
        {
          exp: '(?i)^ganglion interneuron$',
          key: 'ganglion-interneuron',
        },
        {
          exp: '(?i)^giant bipolar cell$',
          key: 'giant-bipolar-cell',
        },
        {
          exp: '(?i)^glandular cell of esophagus$',
          key: 'glandular-cell-of-esophagus',
        },
        {
          exp: '(?i)^glandular epithelial cell$',
          key: 'glandular-epithelial-cell',
        },
        {
          exp: '(?i)^glioblast$',
          key: 'glioblast',
        },
        {
          exp: '(?i)^glomerular capillary endothelial cell$',
          key: 'glomerular-capillary-endothelial-cell',
        },
        {
          exp: '(?i)^glomerular endothelial cell$',
          key: 'glomerular-endothelial-cell',
        },
        {
          exp: '(?i)^glutamatergic neuron$',
          key: 'glutamatergic-neuron',
        },
        {
          exp: '(?i)^glycinergic amacrine cell$',
          key: 'glycinergic-amacrine-cell',
        },
        {
          exp: '(?i)^granular cell of epidermis$',
          key: 'granular-cell-of-epidermis',
        },
        {
          exp: '(?i)^granule cell$',
          key: 'granule-cell',
        },
        {
          exp: '(?i)^granulocyte monocyte progenitor cell$',
          key: 'granulocyte-monocyte-progenitor-cell',
        },
        {
          exp: '(?i)^granulocyte$',
          key: 'granulocyte',
        },
        {
          exp: '(?i)^granulocytopoietic cell$',
          key: 'granulocytopoietic-cell',
        },
        {
          exp: '(?i)^granulosa cell$',
          key: 'granulosa-cell',
        },
        {
          exp: '(?i)^group 2 innate lymphoid cell, human$',
          key: 'group-2-innate-lymphoid-cell-human',
        },
        {
          exp: '(?i)^group 2 innate lymphoid cell$',
          key: 'group-2-innate-lymphoid-cell',
        },
        {
          exp: '(?i)^group 3 innate lymphoid cell, human$',
          key: 'group-3-innate-lymphoid-cell-human',
        },
        {
          exp: '(?i)^group 3 innate lymphoid cell$',
          key: 'group-3-innate-lymphoid-cell',
        },
        {
          exp: '(?i)^gut absorptive cell$',
          key: 'gut-absorptive-cell',
        },
        {
          exp: '(?i)^gut endothelial cell$',
          key: 'gut-endothelial-cell',
        },
        {
          exp: '(?i)^hair follicular keratinocyte$',
          key: 'hair-follicular-keratinocyte',
        },
        {
          exp: '(?i)^helper T cell$',
          key: 'helper-t-cell',
        },
        {
          exp: '(?i)^hematopoietic cell$',
          key: 'hematopoietic-cell',
        },
        {
          exp: '(?i)^hematopoietic multipotent progenitor cell$',
          key: 'hematopoietic-multipotent-progenitor-cell',
        },
        {
          exp: '(?i)^hematopoietic precursor cell$',
          key: 'hematopoietic-precursor-cell',
        },
        {
          exp: '(?i)^hematopoietic stem cell$',
          key: 'hematopoietic-stem-cell',
        },
        {
          exp: '(?i)^hepatic pit cell$',
          key: 'hepatic-pit-cell',
        },
        {
          exp: '(?i)^hepatic stellate cell$',
          key: 'hepatic-stellate-cell',
        },
        {
          exp: '(?i)^hepatoblast$',
          key: 'hepatoblast',
        },
        {
          exp: '(?i)^hippocampal neuron$',
          key: 'hippocampal-neuron',
        },
        {
          exp: '(?i)^hypendymal cell$',
          key: 'hypendymal-cell',
        },
        {
          exp: '(?i)^ileal goblet cell$',
          key: 'ileal-goblet-cell',
        },
        {
          exp: '(?i)^immature B cell$',
          key: 'immature-b-cell',
        },
        {
          exp: '(?i)^immature NK T cell$',
          key: 'immature-nk-t-cell',
        },
        {
          exp: '(?i)^immature Schwann cell$',
          key: 'immature-schwann-cell',
        },
        {
          exp: '(?i)^immature T cell$',
          key: 'immature-t-cell',
        },
        {
          exp: '(?i)^immature alpha-beta T cell$',
          key: 'immature-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^immature astrocyte$',
          key: 'immature-astrocyte',
        },
        {
          exp: '(?i)^immature innate lymphoid cell$',
          key: 'immature-innate-lymphoid-cell',
        },
        {
          exp: '(?i)^immature natural killer cell$',
          key: 'immature-natural-killer-cell',
        },
        {
          exp: '(?i)^immature neutrophil$',
          key: 'immature-neutrophil',
        },
        {
          exp: '(?i)^indirect pathway medium spiny neuron$',
          key: 'indirect-pathway-medium-spiny-neuron',
        },
        {
          exp: '(?i)^direct pathway medium spiny neuron$',
          key: 'direct-pathway-medium-spiny-neuron',
        },
        {
          exp: '(?i)^inflammatory cell$',
          key: 'inflammatory-cell',
        },
        {
          exp: '(?i)^inflammatory macrophage$',
          key: 'inflammatory-macrophage',
        },
        {
          exp: '(?i)^inhibitory interneuron$',
          key: 'inhibitory-interneuron',
        },
        {
          exp: '(?i)^inhibitory motor neuron$',
          key: 'inhibitory-motor-neuron',
        },
        {
          exp: '(?i)^innate lymphoid cell$',
          key: 'innate-lymphoid-cell',
        },
        {
          exp: '(?i)^intermediate monocyte$',
          key: 'intermediate-monocyte',
        },
        {
          exp: '(?i)^interstitial cell of Cajal$',
          key: 'interstitial-cell-of-cajal',
        },
        {
          exp: '(?i)^interstitial cell of ovary$',
          key: 'interstitial-cell-of-ovary',
        },
        {
          exp: '(?i)^intestinal crypt stem cell of colon$',
          key: 'intestinal-crypt-stem-cell-of-colon',
        },
        {
          exp: '(?i)^intestinal crypt stem cell of large intestine$',
          key: 'intestinal-crypt-stem-cell-of-large-intestine',
        },
        {
          exp: '(?i)^intestinal crypt stem cell of small intestine$',
          key: 'intestinal-crypt-stem-cell-of-small-intestine',
        },
        {
          exp: '(?i)^intestinal crypt stem cell$',
          key: 'intestinal-crypt-stem-cell',
        },
        {
          exp: '(?i)^intestinal enteroendocrine cell$',
          key: 'intestinal-enteroendocrine-cell',
        },
        {
          exp: '(?i)^intestinal epithelial cell$',
          key: 'intestinal-epithelial-cell',
        },
        {
          exp: '(?i)^intestinal tuft cell$',
          key: 'intestinal-tuft-cell',
        },
        {
          exp: '(?i)^intraepithelial lymphocyte$',
          key: 'intraepithelial-lymphocyte',
        },
        {
          exp: '(?i)^intrahepatic cholangiocyte$',
          key: 'intrahepatic-cholangiocyte',
        },
        {
          exp: '(?i)^cholangiocyte$',
          key: 'cholangiocyte',
        },
        {
          exp: '(?i)^invaginating midget bipolar cell$',
          key: 'invaginating-midget-bipolar-cell',
        },
        {
          exp: '(?i)^keratinocyte stem cell$',
          key: 'keratinocyte-stem-cell',
        },
        {
          exp: '(?i)^keratocyte$',
          key: 'keratocyte',
        },
        {
          exp: '(?i)^kidney afferent arteriole endothelial cell$',
          key: 'kidney-afferent-arteriole-endothelial-cell',
        },
        {
          exp: '(?i)^kidney capillary endothelial cell$',
          key: 'kidney-capillary-endothelial-cell',
        },
        {
          exp: '(?i)^kidney cell$',
          key: 'kidney-cell',
        },
        {
          exp: '(?i)^kidney collecting duct cell$',
          key: 'kidney-collecting-duct-cell',
        },
        {
          exp: '(?i)^kidney collecting duct epithelial cell$',
          key: 'kidney-collecting-duct-epithelial-cell',
        },
        {
          exp: '(?i)^duct epithelial cell$',
          key: 'duct-epithelial-cell',
        },
        {
          exp: '(?i)^kidney collecting duct intercalated cell$',
          key: 'kidney-collecting-duct-intercalated-cell',
        },
        {
          exp: '(?i)^kidney collecting duct principal cell$',
          key: 'kidney-collecting-duct-principal-cell',
        },
        {
          exp: '(?i)^kidney connecting tubule epithelial cell$',
          key: 'kidney-connecting-tubule-epithelial-cell',
        },
        {
          exp: '(?i)^kidney connecting tubule principal cell$',
          key: 'kidney-connecting-tubule-principal-cell',
        },
        {
          exp: '(?i)^kidney cortex artery cell$',
          key: 'kidney-cortex-artery-cell',
        },
        {
          exp: '(?i)^kidney cortex tubule cell$',
          key: 'kidney-cortex-tubule-cell',
        },
        {
          exp: '(?i)^kidney distal convoluted tubule epithelial cell$',
          key: 'kidney-distal-convoluted-tubule-epithelial-cell',
        },
        {
          exp: '(?i)^kidney efferent arteriole endothelial cell$',
          key: 'kidney-efferent-arteriole-endothelial-cell',
        },
        {
          exp: '(?i)^kidney epithelial cell$',
          key: 'kidney-epithelial-cell',
        },
        {
          exp: '(?i)^kidney glomerular epithelial cell$',
          key: 'kidney-glomerular-epithelial-cell',
        },
        {
          exp: '(?i)^kidney interstitial alternatively activated macrophage$',
          key: 'kidney-interstitial-alternatively-activated-macrophage',
        },
        {
          exp: '(?i)^alternatively activated macrophage$',
          key: 'alternatively-activated-macrophage',
        },
        {
          exp: '(?i)^kidney interstitial cell$',
          key: 'kidney-interstitial-cell',
        },
        {
          exp: '(?i)^kidney interstitial fibroblast$',
          key: 'kidney-interstitial-fibroblast',
        },
        {
          exp: '(?i)^kidney loop of Henle ascending limb epithelial cell$',
          key: 'kidney-loop-of-henle-ascending-limb-epithelial-cell',
        },
        {
          exp: '(?i)^kidney loop of Henle cortical thick ascending limb epithelial cell$',
          key: 'kidney-loop-of-henle-cortical-thick-ascending-limb-epithelial-cell',
        },
        {
          exp: '(?i)^kidney loop of Henle epithelial cell$',
          key: 'kidney-loop-of-henle-epithelial-cell',
        },
        {
          exp: '(?i)^kidney loop of Henle medullary thick ascending limb epithelial cell$',
          key: 'kidney-loop-of-henle-medullary-thick-ascending-limb-epithelial-cell',
        },
        {
          exp: '(?i)^kidney loop of Henle thick ascending limb epithelial cell$',
          key: 'kidney-loop-of-henle-thick-ascending-limb-epithelial-cell',
        },
        {
          exp: '(?i)^kidney loop of Henle thin ascending limb epithelial cell$',
          key: 'kidney-loop-of-henle-thin-ascending-limb-epithelial-cell',
        },
        {
          exp: '(?i)^kidney loop of Henle thin descending limb epithelial cell$',
          key: 'kidney-loop-of-henle-thin-descending-limb-epithelial-cell',
        },
        {
          exp: '(?i)^kidney pelvis urothelial cell$',
          key: 'kidney-pelvis-urothelial-cell',
        },
        {
          exp: '(?i)^kidney proximal convoluted tubule epithelial cell$',
          key: 'kidney-proximal-convoluted-tubule-epithelial-cell',
        },
        {
          exp: '(?i)^kidney proximal straight tubule epithelial cell$',
          key: 'kidney-proximal-straight-tubule-epithelial-cell',
        },
        {
          exp: '(?i)^kidney resident macrophage$',
          key: 'kidney-resident-macrophage',
        },
        {
          exp: '(?i)^lactocyte$',
          key: 'lactocyte',
        },
        {
          exp: '(?i)^lamp5 GABAergic cortical interneuron$',
          key: 'lamp5-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^large intestine goblet cell$',
          key: 'large-intestine-goblet-cell',
        },
        {
          exp: '(?i)^large pre-B-II cell$',
          key: 'large-pre-b-ii-cell',
        },
        {
          exp: '(?i)^late pro-B cell$',
          key: 'late-pro-b-cell',
        },
        {
          exp: '(?i)^late promyelocyte$',
          key: 'late-promyelocyte',
        },
        {
          exp: '(?i)^lens epithelial cell$',
          key: 'lens-epithelial-cell',
        },
        {
          exp: '(?i)^lens fiber cell$',
          key: 'lens-fiber-cell',
        },
        {
          exp: '(?i)^liver dendritic cell$',
          key: 'liver-dendritic-cell',
        },
        {
          exp: '(?i)^luminal adaptive secretory precursor cell of mammary gland$',
          key: 'luminal-adaptive-secretory-precursor-cell-of-mammary-gland',
        },
        {
          exp: '(?i)^luminal cell of prostate epithelium$',
          key: 'luminal-cell-of-prostate-epithelium',
        },
        {
          exp: '(?i)^luminal epithelial cell of mammary gland$',
          key: 'luminal-epithelial-cell-of-mammary-gland',
        },
        {
          exp: '(?i)^luminal hormone-sensing cell of mammary gland$',
          key: 'luminal-hormone-sensing-cell-of-mammary-gland',
        },
        {
          exp: '(?i)^lung ciliated cell$',
          key: 'lung-ciliated-cell',
        },
        {
          exp: '(?i)^ciliated cell$',
          key: 'ciliated-cell',
        },
        {
          exp: '(?i)^lung goblet cell$',
          key: 'lung-goblet-cell',
        },
        {
          exp: '(?i)^lung interstitial macrophage$',
          key: 'lung-interstitial-macrophage',
        },
        {
          exp: '(?i)^lung macrophage$',
          key: 'lung-macrophage',
        },
        {
          exp: '(?i)^lung microvascular endothelial cell$',
          key: 'lung-microvascular-endothelial-cell',
        },
        {
          exp: '(?i)^lung neuroendocrine cell$',
          key: 'lung-neuroendocrine-cell',
        },
        {
          exp: '(?i)^lung perichondrial fibroblast$',
          key: 'lung-perichondrial-fibroblast',
        },
        {
          exp: '(?i)^lung pericyte$',
          key: 'lung-pericyte',
        },
        {
          exp: '(?i)^lung resident memory CD4-positive, alpha-beta T cell$',
          key: 'lung-resident-memory-cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^lung resident memory CD8-positive, alpha-beta T cell$',
          key: 'lung-resident-memory-cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^lung secretory cell$',
          key: 'lung-secretory-cell',
        },
        {
          exp: '(?i)^lymphocyte of B lineage$',
          key: 'lymphocyte-of-b-lineage',
        },
        {
          exp: '(?i)^lymphocyte$',
          key: 'lymphocyte',
        },
        {
          exp: '(?i)^lymphoid lineage restricted progenitor cell$',
          key: 'lymphoid-lineage-restricted-progenitor-cell',
        },
        {
          exp: '(?i)^macroglial cell$',
          key: 'macroglial-cell',
        },
        {
          exp: '(?i)^macrophage dendritic cell progenitor$',
          key: 'macrophage-dendritic-cell-progenitor',
        },
        {
          exp: '(?i)^macula densa epithelial cell$',
          key: 'macula-densa-epithelial-cell',
        },
        {
          exp: '(?i)^male germ cell$',
          key: 'male-germ-cell',
        },
        {
          exp: '(?i)^malignant cell$',
          key: 'malignant-cell',
        },
        {
          exp: '(?i)^mammary gland epithelial cell$',
          key: 'mammary-gland-epithelial-cell',
        },
        {
          exp: '(?i)^mast cell$',
          key: 'mast-cell',
        },
        {
          exp: '(?i)^mature B cell$',
          key: 'mature-b-cell',
        },
        {
          exp: '(?i)^mature NK T cell$',
          key: 'mature-nk-t-cell',
        },
        {
          exp: '(?i)^mature T cell$',
          key: 'mature-t-cell',
        },
        {
          exp: '(?i)^mature alpha-beta T cell$',
          key: 'mature-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^mature astrocyte$',
          key: 'mature-astrocyte',
        },
        {
          exp: '(?i)^mature conventional dendritic cell$',
          key: 'mature-conventional-dendritic-cell',
        },
        {
          exp: '(?i)^mature gamma-delta T cell$',
          key: 'mature-gamma-delta-t-cell',
        },
        {
          exp: '(?i)^gamma-delta T cell$',
          key: 'gamma-delta-t-cell',
        },
        {
          exp: '(?i)^mature microglial cell$',
          key: 'mature-microglial-cell',
        },
        {
          exp: '(?i)^medium spiny neuron$',
          key: 'medium-spiny-neuron',
        },
        {
          exp: '(?i)^medullary thymic epithelial cell$',
          key: 'medullary-thymic-epithelial-cell',
        },
        {
          exp: '(?i)^megakaryocyte progenitor cell$',
          key: 'megakaryocyte-progenitor-cell',
        },
        {
          exp: '(?i)^megakaryocyte-erythroid progenitor cell$',
          key: 'megakaryocyte-erythroid-progenitor-cell',
        },
        {
          exp: '(?i)^erythroid progenitor cell$',
          key: 'erythroid-progenitor-cell',
        },
        {
          exp: '(?i)^megakaryocyte$',
          key: 'megakaryocyte',
        },
        {
          exp: '(?i)^meis2 expressing cortical GABAergic cell$',
          key: 'meis2-expressing-cortical-gabaergic-cell',
        },
        {
          exp: '(?i)^melanocyte of skin$',
          key: 'melanocyte-of-skin',
        },
        {
          exp: '(?i)^melanocyte$',
          key: 'melanocyte',
        },
        {
          exp: '(?i)^memory T cell$',
          key: 'memory-t-cell',
        },
        {
          exp: '(?i)^memory regulatory T cell$',
          key: 'memory-regulatory-t-cell',
        },
        {
          exp: '(?i)^meningeal macrophage$',
          key: 'meningeal-macrophage',
        },
        {
          exp: '(?i)^mesangial cell$',
          key: 'mesangial-cell',
        },
        {
          exp: '(?i)^mesenchymal cell$',
          key: 'mesenchymal-cell',
        },
        {
          exp: '(?i)^mesenchymal lymphangioblast$',
          key: 'mesenchymal-lymphangioblast',
        },
        {
          exp: '(?i)^mesenchymal stem cell of adipose tissue$',
          key: 'mesenchymal-stem-cell-of-adipose-tissue',
        },
        {
          exp: '(?i)^mesenchymal stem cell$',
          key: 'mesenchymal-stem-cell',
        },
        {
          exp: '(?i)^mesodermal cell$',
          key: 'mesodermal-cell',
        },
        {
          exp: '(?i)^mesothelial cell of epicardium$',
          key: 'mesothelial-cell-of-epicardium',
        },
        {
          exp: '(?i)^mesothelial cell of pleura$',
          key: 'mesothelial-cell-of-pleura',
        },
        {
          exp: '(?i)^mesothelial cell of visceral pleura$',
          key: 'mesothelial-cell-of-visceral-pleura',
        },
        {
          exp: '(?i)^mesothelial cell$',
          key: 'mesothelial-cell',
        },
        {
          exp: '(?i)^mesothelial fibroblast$',
          key: 'mesothelial-fibroblast',
        },
        {
          exp: '(?i)^metallothionein-positive alveolar macrophage$',
          key: 'metallothionein-positive-alveolar-macrophage',
        },
        {
          exp: '(?i)^alveolar macrophage$',
          key: 'alveolar-macrophage',
        },
        {
          exp: '(?i)^microfold cell of epithelium of small intestine$',
          key: 'microfold-cell-of-epithelium-of-small-intestine',
        },
        {
          exp: '(?i)^microglial cell$',
          key: 'microglial-cell',
        },
        {
          exp: '(?i)^midget ganglion cell of retina$',
          key: 'midget-ganglion-cell-of-retina',
        },
        {
          exp: '(?i)^midzonal region hepatocyte$',
          key: 'midzonal-region-hepatocyte',
        },
        {
          exp: '(?i)^migratory enteric neural crest cell$',
          key: 'migratory-enteric-neural-crest-cell',
        },
        {
          exp: '(?i)^monocyte-derived dendritic cell$',
          key: 'monocyte-derived-dendritic-cell',
        },
        {
          exp: '(?i)^mononuclear phagocyte$',
          key: 'mononuclear-phagocyte',
        },
        {
          exp: '(?i)^mucosal invariant T cell$',
          key: 'mucosal-invariant-t-cell',
        },
        {
          exp: '(?i)^mucous neck cell$',
          key: 'mucous-neck-cell',
        },
        {
          exp: '(?i)^mucus secreting cell$',
          key: 'mucus-secreting-cell',
        },
        {
          exp: '(?i)^multi-ciliated epithelial cell$',
          key: 'multi-ciliated-epithelial-cell',
        },
        {
          exp: '(?i)^ciliated epithelial cell$',
          key: 'ciliated-epithelial-cell',
        },
        {
          exp: '(?i)^mural cell$',
          key: 'mural-cell',
        },
        {
          exp: '(?i)^muscle precursor cell$',
          key: 'muscle-precursor-cell',
        },
        {
          exp: '(?i)^myeloid cell$',
          key: 'myeloid-cell',
        },
        {
          exp: '(?i)^myeloid dendritic cell, human$',
          key: 'myeloid-dendritic-cell-human',
        },
        {
          exp: '(?i)^myeloid dendritic cell$',
          key: 'myeloid-dendritic-cell',
        },
        {
          exp: '(?i)^myeloid leukocyte$',
          key: 'myeloid-leukocyte',
        },
        {
          exp: '(?i)^leukocyte$',
          key: 'leukocyte',
        },
        {
          exp: '(?i)^myeloid lineage restricted progenitor cell$',
          key: 'myeloid-lineage-restricted-progenitor-cell',
        },
        {
          exp: '(?i)^myoepithelial cell of mammary gland$',
          key: 'myoepithelial-cell-of-mammary-gland',
        },
        {
          exp: '(?i)^myoepithelial cell$',
          key: 'myoepithelial-cell',
        },
        {
          exp: '(?i)^myofibroblast cell$',
          key: 'myofibroblast-cell',
        },
        {
          exp: '(?i)^myometrial cell$',
          key: 'myometrial-cell',
        },
        {
          exp: '(?i)^naive B cell$',
          key: 'naive-b-cell',
        },
        {
          exp: '(?i)^naive T cell$',
          key: 'naive-t-cell',
        },
        {
          exp: '(?i)^naive regulatory T cell$',
          key: 'naive-regulatory-t-cell',
        },
        {
          exp: '(?i)^naive thymus-derived CD4-positive, alpha-beta T cell$',
          key: 'naive-thymus-derived-cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^CD4-positive, alpha-beta T cell$',
          key: 'cd4-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^naive thymus-derived CD8-positive, alpha-beta T cell$',
          key: 'naive-thymus-derived-cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^CD8-positive, alpha-beta T cell$',
          key: 'cd8-positive-alpha-beta-t-cell',
        },
        {
          exp: '(?i)^alpha-beta T cell$',
          key: 'alpha-beta-t-cell',
        },
        {
          exp: '(?i)^nasal mucosa goblet cell$',
          key: 'nasal-mucosa-goblet-cell',
        },
        {
          exp: '(?i)^natural T-regulatory cell$',
          key: 'natural-t-regulatory-cell',
        },
        {
          exp: '(?i)^near-projecting glutamatergic cortical neuron$',
          key: 'near-projecting-glutamatergic-cortical-neuron',
        },
        {
          exp: '(?i)^neoplastic cell$',
          key: 'neoplastic-cell',
        },
        {
          exp: '(?i)^nephron tubule epithelial cell$',
          key: 'nephron-tubule-epithelial-cell',
        },
        {
          exp: '(?i)^neural cell$',
          key: 'neural-cell',
        },
        {
          exp: '(?i)^neural crest cell$',
          key: 'neural-crest-cell',
        },
        {
          exp: '(?i)^neural progenitor cell$',
          key: 'neural-progenitor-cell',
        },
        {
          exp: '(?i)^neural stem cell$',
          key: 'neural-stem-cell',
        },
        {
          exp: '(?i)^neuroblast \\(sensu Vertebrata\\)$',
          key: 'neuroblast-sensu-vertebrata',
        },
        {
          exp: '(?i)^neuroendocrine cell$',
          key: 'neuroendocrine-cell',
        },
        {
          exp: '(?i)^neuroepithelial stem cell$',
          key: 'neuroepithelial-stem-cell',
        },
        {
          exp: '(?i)^neuron associated cell \\(sensu Vertebrata\\)$',
          key: 'neuron-associated-cell-sensu-vertebrata',
        },
        {
          exp: '(?i)^neuron associated cell$',
          key: 'neuron-associated-cell',
        },
        {
          exp: '(?i)^neuronal brush cell$',
          key: 'neuronal-brush-cell',
        },
        {
          exp: '(?i)^neuronal receptor cell$',
          key: 'neuronal-receptor-cell',
        },
        {
          exp: '(?i)^neutrophil progenitor cell$',
          key: 'neutrophil-progenitor-cell',
        },
        {
          exp: '(?i)^neutrophil$',
          key: 'neutrophil',
        },
        {
          exp: '(?i)^non-classical monocyte$',
          key: 'non-classical-monocyte',
        },
        {
          exp: '(?i)^classical monocyte$',
          key: 'classical-monocyte',
        },
        {
          exp: '(?i)^non-myelinating Schwann cell$',
          key: 'non-myelinating-schwann-cell',
        },
        {
          exp: '(?i)^myelinating Schwann cell$',
          key: 'myelinating-schwann-cell',
        },
        {
          exp: '(?i)^Schwann cell$',
          key: 'schwann-cell',
        },
        {
          exp: '(?i)^non-pigmented ciliary epithelial cell$',
          key: 'non-pigmented-ciliary-epithelial-cell',
        },
        {
          exp: '(?i)^non-terminally differentiated cell$',
          key: 'non-terminally-differentiated-cell',
        },
        {
          exp: '(?i)^noradrenergic cell$',
          key: 'noradrenergic-cell',
        },
        {
          exp: '(?i)^oligodendrocyte precursor cell$',
          key: 'oligodendrocyte-precursor-cell',
        },
        {
          exp: '(?i)^oligodendrocyte$',
          key: 'oligodendrocyte',
        },
        {
          exp: '(?i)^oocyte$',
          key: 'oocyte',
        },
        {
          exp: '(?i)^oogonial cell$',
          key: 'oogonial-cell',
        },
        {
          exp: '(?i)^osteoclast$',
          key: 'osteoclast',
        },
        {
          exp: '(?i)^ovarian surface epithelial cell$',
          key: 'ovarian-surface-epithelial-cell',
        },
        {
          exp: '(?i)^pancreatic A cell$',
          key: 'pancreatic-a-cell',
        },
        {
          exp: '(?i)^pancreatic D cell$',
          key: 'pancreatic-d-cell',
        },
        {
          exp: '(?i)^pancreatic PP cell$',
          key: 'pancreatic-pp-cell',
        },
        {
          exp: '(?i)^PP cell$',
          key: 'pp-cell',
        },
        {
          exp: '(?i)^pancreatic acinar cell$',
          key: 'pancreatic-acinar-cell',
        },
        {
          exp: '(?i)^acinar cell$',
          key: 'acinar-cell',
        },
        {
          exp: '(?i)^pancreatic ductal cell$',
          key: 'pancreatic-ductal-cell',
        },
        {
          exp: '(?i)^pancreatic endocrine cell$',
          key: 'pancreatic-endocrine-cell',
        },
        {
          exp: '(?i)^pancreatic epsilon cell$',
          key: 'pancreatic-epsilon-cell',
        },
        {
          exp: '(?i)^pancreatic stellate cell$',
          key: 'pancreatic-stellate-cell',
        },
        {
          exp: '(?i)^paneth cell of colon$',
          key: 'paneth-cell-of-colon',
        },
        {
          exp: '(?i)^paneth cell of epithelium of small intestine$',
          key: 'paneth-cell-of-epithelium-of-small-intestine',
        },
        {
          exp: '(?i)^paneth cell$',
          key: 'paneth-cell',
        },
        {
          exp: '(?i)^papillary tips cell$',
          key: 'papillary-tips-cell',
        },
        {
          exp: '(?i)^parasol ganglion cell of retina$',
          key: 'parasol-ganglion-cell-of-retina',
        },
        {
          exp: '(?i)^parietal cell$',
          key: 'parietal-cell',
        },
        {
          exp: '(?i)^parietal epithelial cell$',
          key: 'parietal-epithelial-cell',
        },
        {
          exp: '(?i)^peptic cell$',
          key: 'peptic-cell',
        },
        {
          exp: '(?i)^peripheral blood mononuclear cell$',
          key: 'peripheral-blood-mononuclear-cell',
        },
        {
          exp: '(?i)^mononuclear cell$',
          key: 'mononuclear-cell',
        },
        {
          exp: '(?i)^peripheral nervous system neuron$',
          key: 'peripheral-nervous-system-neuron',
        },
        {
          exp: '(?i)^periportal region hepatocyte$',
          key: 'periportal-region-hepatocyte',
        },
        {
          exp: '(?i)^hepatocyte$',
          key: 'hepatocyte',
        },
        {
          exp: '(?i)^peritubular capillary endothelial cell$',
          key: 'peritubular-capillary-endothelial-cell',
        },
        {
          exp: '(?i)^capillary endothelial cell$',
          key: 'capillary-endothelial-cell',
        },
        {
          exp: '(?i)^perivascular cell$',
          key: 'perivascular-cell',
        },
        {
          exp: '(?i)^phagocyte$',
          key: 'phagocyte',
        },
        {
          exp: '(?i)^photoreceptor cell$',
          key: 'photoreceptor-cell',
        },
        {
          exp: '(?i)^pigmented ciliary epithelial cell$',
          key: 'pigmented-ciliary-epithelial-cell',
        },
        {
          exp: '(?i)^pigmented epithelial cell$',
          key: 'pigmented-epithelial-cell',
        },
        {
          exp: '(?i)^pituitary gland cell$',
          key: 'pituitary-gland-cell',
        },
        {
          exp: '(?i)^placental villous trophoblast$',
          key: 'placental-villous-trophoblast',
        },
        {
          exp: '(?i)^plasma cell$',
          key: 'plasma-cell',
        },
        {
          exp: '(?i)^plasmablast$',
          key: 'plasmablast',
        },
        {
          exp: '(?i)^plasmacytoid dendritic cell, human$',
          key: 'plasmacytoid-dendritic-cell-human',
        },
        {
          exp: '(?i)^dendritic cell, human$',
          key: 'dendritic-cell-human',
        },
        {
          exp: '(?i)^plasmacytoid dendritic cell$',
          key: 'plasmacytoid-dendritic-cell',
        },
        {
          exp: '(?i)^plasmatocyte$',
          key: 'plasmatocyte',
        },
        {
          exp: '(?i)^platelet$',
          key: 'platelet',
        },
        {
          exp: '(?i)^podocyte$',
          key: 'podocyte',
        },
        {
          exp: '(?i)^pre-B-I cell$',
          key: 'pre-b-i-cell',
        },
        {
          exp: '(?i)^pre-conventional dendritic cell$',
          key: 'pre-conventional-dendritic-cell',
        },
        {
          exp: '(?i)^conventional dendritic cell$',
          key: 'conventional-dendritic-cell',
        },
        {
          exp: '(?i)^dendritic cell$',
          key: 'dendritic-cell',
        },
        {
          exp: '(?i)^pre-natural killer cell$',
          key: 'pre-natural-killer-cell',
        },
        {
          exp: '(?i)^natural killer cell$',
          key: 'natural-killer-cell',
        },
        {
          exp: '(?i)^preadipocyte$',
          key: 'preadipocyte',
        },
        {
          exp: '(?i)^precursor B cell$',
          key: 'precursor-b-cell',
        },
        {
          exp: '(?i)^precursor cell$',
          key: 'precursor-cell',
        },
        {
          exp: '(?i)^preosteoblast$',
          key: 'preosteoblast',
        },
        {
          exp: '(?i)^osteoblast$',
          key: 'osteoblast',
        },
        {
          exp: '(?i)^prickle cell$',
          key: 'prickle-cell',
        },
        {
          exp: '(?i)^primary cultured cell$',
          key: 'primary-cultured-cell',
        },
        {
          exp: '(?i)^primary sensory neuron \\(sensu Teleostei\\)$',
          key: 'primary-sensory-neuron-sensu-teleostei',
        },
        {
          exp: '(?i)^primitive red blood cell$',
          key: 'primitive-red-blood-cell',
        },
        {
          exp: '(?i)^blood cell$',
          key: 'blood-cell',
        },
        {
          exp: '(?i)^primordial germ cell$',
          key: 'primordial-germ-cell',
        },
        {
          exp: '(?i)^germ cell$',
          key: 'germ-cell',
        },
        {
          exp: '(?i)^pro-B cell$',
          key: 'pro-b-cell',
        },
        {
          exp: '(?i)^pro-T cell$',
          key: 'pro-t-cell',
        },
        {
          exp: '(?i)^proerythroblast$',
          key: 'proerythroblast',
        },
        {
          exp: '(?i)^erythroblast$',
          key: 'erythroblast',
        },
        {
          exp: '(?i)^professional antigen presenting cell$',
          key: 'professional-antigen-presenting-cell',
        },
        {
          exp: '(?i)^progenitor cell of endocrine pancreas$',
          key: 'progenitor-cell-of-endocrine-pancreas',
        },
        {
          exp: '(?i)^progenitor cell of mammary luminal epithelium$',
          key: 'progenitor-cell-of-mammary-luminal-epithelium',
        },
        {
          exp: '(?i)^progenitor cell$',
          key: 'progenitor-cell',
        },
        {
          exp: '(?i)^promonocyte$',
          key: 'promonocyte',
        },
        {
          exp: '(?i)^monocyte$',
          key: 'monocyte',
        },
        {
          exp: '(?i)^promyelocyte$',
          key: 'promyelocyte',
        },
        {
          exp: '(?i)^myelocyte$',
          key: 'myelocyte',
        },
        {
          exp: '(?i)^prostate gland microvascular endothelial cell$',
          key: 'prostate-gland-microvascular-endothelial-cell',
        },
        {
          exp: '(?i)^microvascular endothelial cell$',
          key: 'microvascular-endothelial-cell',
        },
        {
          exp: '(?i)^prostate stromal cell$',
          key: 'prostate-stromal-cell',
        },
        {
          exp: '(?i)^pulmonary artery endothelial cell$',
          key: 'pulmonary-artery-endothelial-cell',
        },
        {
          exp: '(?i)^pulmonary interstitial fibroblast$',
          key: 'pulmonary-interstitial-fibroblast',
        },
        {
          exp: '(?i)^pulmonary ionocyte$',
          key: 'pulmonary-ionocyte',
        },
        {
          exp: '(?i)^ionocyte$',
          key: 'ionocyte',
        },
        {
          exp: '(?i)^pvalb GABAergic cortical interneuron$',
          key: 'pvalb-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^pyramidal neuron$',
          key: 'pyramidal-neuron',
        },
        {
          exp: '(?i)^radial glial cell$',
          key: 'radial-glial-cell',
        },
        {
          exp: '(?i)^glial cell$',
          key: 'glial-cell',
        },
        {
          exp: '(?i)^regular atrial cardiac myocyte$',
          key: 'regular-atrial-cardiac-myocyte',
        },
        {
          exp: '(?i)^regular ventricular cardiac myocyte$',
          key: 'regular-ventricular-cardiac-myocyte',
        },
        {
          exp: '(?i)^regulatory T cell$',
          key: 'regulatory-t-cell',
        },
        {
          exp: '(?i)^renal alpha-intercalated cell$',
          key: 'renal-alpha-intercalated-cell',
        },
        {
          exp: '(?i)^renal beta-intercalated cell$',
          key: 'renal-beta-intercalated-cell',
        },
        {
          exp: '(?i)^renal intercalated cell$',
          key: 'renal-intercalated-cell',
        },
        {
          exp: '(?i)^renal interstitial pericyte$',
          key: 'renal-interstitial-pericyte',
        },
        {
          exp: '(?i)^pericyte$',
          key: 'pericyte',
        },
        {
          exp: '(?i)^renal principal cell$',
          key: 'renal-principal-cell',
        },
        {
          exp: '(?i)^respiratory basal cell$',
          key: 'respiratory-basal-cell',
        },
        {
          exp: '(?i)^respiratory epithelial cell$',
          key: 'respiratory-epithelial-cell',
        },
        {
          exp: '(?i)^respiratory goblet cell$',
          key: 'respiratory-goblet-cell',
        },
        {
          exp: '(?i)^respiratory hillock cell$',
          key: 'respiratory-hillock-cell',
        },
        {
          exp: '(?i)^respiratory suprabasal cell$',
          key: 'respiratory-suprabasal-cell',
        },
        {
          exp: '(?i)^basal cell$',
          key: 'basal-cell',
        },
        {
          exp: '(?i)^reticular cell$',
          key: 'reticular-cell',
        },
        {
          exp: '(?i)^reticulocyte$',
          key: 'reticulocyte',
        },
        {
          exp: '(?i)^retina horizontal cell$',
          key: 'retina-horizontal-cell',
        },
        {
          exp: '(?i)^retinal astrocyte$',
          key: 'retinal-astrocyte',
        },
        {
          exp: '(?i)^astrocyte$',
          key: 'astrocyte',
        },
        {
          exp: '(?i)^retinal bipolar neuron$',
          key: 'retinal-bipolar-neuron',
        },
        {
          exp: '(?i)^bipolar neuron$',
          key: 'bipolar-neuron',
        },
        {
          exp: '(?i)^retinal blood vessel endothelial cell$',
          key: 'retinal-blood-vessel-endothelial-cell',
        },
        {
          exp: '(?i)^blood vessel endothelial cell$',
          key: 'blood-vessel-endothelial-cell',
        },
        {
          exp: '(?i)^retinal cone cell$',
          key: 'retinal-cone-cell',
        },
        {
          exp: '(?i)^retinal ganglion cell$',
          key: 'retinal-ganglion-cell',
        },
        {
          exp: '(?i)^retinal pigment epithelial cell$',
          key: 'retinal-pigment-epithelial-cell',
        },
        {
          exp: '(?i)^retinal rod cell$',
          key: 'retinal-rod-cell',
        },
        {
          exp: '(?i)^rod bipolar cell$',
          key: 'rod-bipolar-cell',
        },
        {
          exp: '(?i)^salivary gland cell$',
          key: 'salivary-gland-cell',
        },
        {
          exp: '(?i)^sebaceous gland cell$',
          key: 'sebaceous-gland-cell',
        },
        {
          exp: '(?i)^sebum secreting cell$',
          key: 'sebum-secreting-cell',
        },
        {
          exp: '(?i)^secondary lens fiber$',
          key: 'secondary-lens-fiber',
        },
        {
          exp: '(?i)^secretory cell$',
          key: 'secretory-cell',
        },
        {
          exp: '(?i)^seminal vesicle glandular cell$',
          key: 'seminal-vesicle-glandular-cell',
        },
        {
          exp: '(?i)^sensory neuron$',
          key: 'sensory-neuron',
        },
        {
          exp: '(?i)^serous cell of epithelium of bronchus$',
          key: 'serous-cell-of-epithelium-of-bronchus',
        },
        {
          exp: '(?i)^serous cell of epithelium of trachea$',
          key: 'serous-cell-of-epithelium-of-trachea',
        },
        {
          exp: '(?i)^serous secreting cell of bronchus submucosal gland$',
          key: 'serous-secreting-cell-of-bronchus-submucosal-gland',
        },
        {
          exp: '(?i)^serous secreting cell$',
          key: 'serous-secreting-cell',
        },
        {
          exp: '(?i)^skeletal muscle fiber$',
          key: 'skeletal-muscle-fiber',
        },
        {
          exp: '(?i)^skeletal muscle fibroblast$',
          key: 'skeletal-muscle-fibroblast',
        },
        {
          exp: '(?i)^skeletal muscle satellite cell$',
          key: 'skeletal-muscle-satellite-cell',
        },
        {
          exp: '(?i)^skeletal muscle satellite stem cell$',
          key: 'skeletal-muscle-satellite-stem-cell',
        },
        {
          exp: '(?i)^skin fibroblast$',
          key: 'skin-fibroblast',
        },
        {
          exp: '(?i)^fibroblast$',
          key: 'fibroblast',
        },
        {
          exp: '(?i)^slow muscle cell$',
          key: 'slow-muscle-cell',
        },
        {
          exp: '(?i)^small intestine goblet cell$',
          key: 'small-intestine-goblet-cell',
        },
        {
          exp: '(?i)^intestine goblet cell$',
          key: 'intestine-goblet-cell',
        },
        {
          exp: '(?i)^small pre-B-II cell$',
          key: 'small-pre-b-ii-cell',
        },
        {
          exp: '(?i)^smooth muscle cell of large intestine$',
          key: 'smooth-muscle-cell-of-large-intestine',
        },
        {
          exp: '(?i)^smooth muscle cell of prostate$',
          key: 'smooth-muscle-cell-of-prostate',
        },
        {
          exp: '(?i)^smooth muscle cell of small intestine$',
          key: 'smooth-muscle-cell-of-small-intestine',
        },
        {
          exp: '(?i)^smooth muscle cell of sphincter of pupil$',
          key: 'smooth-muscle-cell-of-sphincter-of-pupil',
        },
        {
          exp: '(?i)^smooth muscle cell of the pulmonary artery$',
          key: 'smooth-muscle-cell-of-the-pulmonary-artery',
        },
        {
          exp: '(?i)^smooth muscle cell of trachea$',
          key: 'smooth-muscle-cell-of-trachea',
        },
        {
          exp: '(?i)^smooth muscle fiber of ileum$',
          key: 'smooth-muscle-fiber-of-ileum',
        },
        {
          exp: '(?i)^smooth muscle myoblast$',
          key: 'smooth-muscle-myoblast',
        },
        {
          exp: '(?i)^myoblast$',
          key: 'myoblast',
        },
        {
          exp: '(?i)^sncg GABAergic cortical interneuron$',
          key: 'sncg-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^sperm$',
          key: 'sperm',
        },
        {
          exp: '(?i)^squamous epithelial cell$',
          key: 'squamous-epithelial-cell',
        },
        {
          exp: '(?i)^sst GABAergic cortical interneuron$',
          key: 'sst-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^sst chodl GABAergic cortical interneuron$',
          key: 'sst-chodl-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^starburst amacrine cell$',
          key: 'starburst-amacrine-cell',
        },
        {
          exp: '(?i)^amacrine cell$',
          key: 'amacrine-cell',
        },
        {
          exp: '(?i)^stellate neuron$',
          key: 'stellate-neuron',
        },
        {
          exp: '(?i)^stem cell of epidermis$',
          key: 'stem-cell-of-epidermis',
        },
        {
          exp: '(?i)^stem cell$',
          key: 'stem-cell',
        },
        {
          exp: '(?i)^stratified epithelial cell$',
          key: 'stratified-epithelial-cell',
        },
        {
          exp: '(?i)^epithelial cell$',
          key: 'epithelial-cell',
        },
        {
          exp: '(?i)^stromal cell of endometrium$',
          key: 'stromal-cell-of-endometrium',
        },
        {
          exp: '(?i)^stromal cell of lamina propria of small intestine$',
          key: 'stromal-cell-of-lamina-propria-of-small-intestine',
        },
        {
          exp: '(?i)^stromal cell of ovary$',
          key: 'stromal-cell-of-ovary',
        },
        {
          exp: '(?i)^stromal cell$',
          key: 'stromal-cell',
        },
        {
          exp: '(?i)^subcutaneous adipocyte$',
          key: 'subcutaneous-adipocyte',
        },
        {
          exp: '(?i)^adipocyte$',
          key: 'adipocyte',
        },
        {
          exp: '(?i)^supporting cell$',
          key: 'supporting-cell',
        },
        {
          exp: '(?i)^suprabasal keratinocyte$',
          key: 'suprabasal-keratinocyte',
        },
        {
          exp: '(?i)^keratinocyte$',
          key: 'keratinocyte',
        },
        {
          exp: '(?i)^surface ectodermal cell$',
          key: 'surface-ectodermal-cell',
        },
        {
          exp: '(?i)^sympathetic neuron$',
          key: 'sympathetic-neuron',
        },
        {
          exp: '(?i)^syncytiotrophoblast cell$',
          key: 'syncytiotrophoblast-cell',
        },
        {
          exp: '(?i)^tanycyte$',
          key: 'tanycyte',
        },
        {
          exp: '(?i)^taste receptor cell$',
          key: 'taste-receptor-cell',
        },
        {
          exp: '(?i)^tendon cell$',
          key: 'tendon-cell',
        },
        {
          exp: '(?i)^theca cell$',
          key: 'theca-cell',
        },
        {
          exp: '(?i)^thymocyte$',
          key: 'thymocyte',
        },
        {
          exp: '(?i)^thyroid follicular cell$',
          key: 'thyroid-follicular-cell',
        },
        {
          exp: '(?i)^tissue-resident macrophage$',
          key: 'tissue-resident-macrophage',
        },
        {
          exp: '(?i)^macrophage$',
          key: 'macrophage',
        },
        {
          exp: '(?i)^tongue muscle cell$',
          key: 'tongue-muscle-cell',
        },
        {
          exp: '(?i)^tonsil germinal center B cell$',
          key: 'tonsil-germinal-center-b-cell',
        },
        {
          exp: '(?i)^germinal center B cell$',
          key: 'germinal-center-b-cell',
        },
        {
          exp: '(?i)^tracheal goblet cell$',
          key: 'tracheal-goblet-cell',
        },
        {
          exp: '(?i)^tracheobronchial goblet cell$',
          key: 'tracheobronchial-goblet-cell',
        },
        {
          exp: '(?i)^bronchial goblet cell$',
          key: 'bronchial-goblet-cell',
        },
        {
          exp: '(?i)^goblet cell$',
          key: 'goblet-cell',
        },
        {
          exp: '(?i)^tracheobronchial serous cell$',
          key: 'tracheobronchial-serous-cell',
        },
        {
          exp: '(?i)^tracheobronchial smooth muscle cell$',
          key: 'tracheobronchial-smooth-muscle-cell',
        },
        {
          exp: '(?i)^bronchial smooth muscle cell$',
          key: 'bronchial-smooth-muscle-cell',
        },
        {
          exp: '(?i)^transit amplifying cell of colon$',
          key: 'transit-amplifying-cell-of-colon',
        },
        {
          exp: '(?i)^transit amplifying cell of small intestine$',
          key: 'transit-amplifying-cell-of-small-intestine',
        },
        {
          exp: '(?i)^transit amplifying cell$',
          key: 'transit-amplifying-cell',
        },
        {
          exp: '(?i)^transitional stage B cell$',
          key: 'transitional-stage-b-cell',
        },
        {
          exp: '(?i)^trophoblast giant cell$',
          key: 'trophoblast-giant-cell',
        },
        {
          exp: '(?i)^tuft cell of colon$',
          key: 'tuft-cell-of-colon',
        },
        {
          exp: '(?i)^type 1 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-1-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 2 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-2-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 3a cone bipolar cell$',
          key: 'type-3a-cone-bipolar-cell',
        },
        {
          exp: '(?i)^type 3b cone bipolar cell$',
          key: 'type-3b-cone-bipolar-cell',
        },
        {
          exp: '(?i)^type 4 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-4-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 5 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-5-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 5a cone bipolar cell$',
          key: 'type-5a-cone-bipolar-cell',
        },
        {
          exp: '(?i)^type 5b cone bipolar cell$',
          key: 'type-5b-cone-bipolar-cell',
        },
        {
          exp: '(?i)^type 6 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-6-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 7 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-7-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 8 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-8-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type 9 cone bipolar cell \\(sensu Mus\\)$',
          key: 'type-9-cone-bipolar-cell-sensu-mus',
        },
        {
          exp: '(?i)^type B pancreatic cell$',
          key: 'type-b-pancreatic-cell',
        },
        {
          exp: '(?i)^type D enteroendocrine cell$',
          key: 'type-d-enteroendocrine-cell',
        },
        {
          exp: '(?i)^type EC enteroendocrine cell$',
          key: 'type-ec-enteroendocrine-cell',
        },
        {
          exp: '(?i)^type G enteroendocrine cell$',
          key: 'type-g-enteroendocrine-cell',
        },
        {
          exp: '(?i)^type I NK T cell$',
          key: 'type-i-nk-t-cell',
        },
        {
          exp: '(?i)^type I enteroendocrine cell$',
          key: 'type-i-enteroendocrine-cell',
        },
        {
          exp: '(?i)^type I pneumocyte$',
          key: 'type-i-pneumocyte',
        },
        {
          exp: '(?i)^type II NK T cell$',
          key: 'type-ii-nk-t-cell',
        },
        {
          exp: '(?i)^T cell$',
          key: 't-cell',
        },
        {
          exp: '(?i)^type II pneumocyte$',
          key: 'type-ii-pneumocyte',
        },
        {
          exp: '(?i)^pneumocyte$',
          key: 'pneumocyte',
        },
        {
          exp: '(?i)^type L enteroendocrine cell$',
          key: 'type-l-enteroendocrine-cell',
        },
        {
          exp: '(?i)^type N enteroendocrine cell$',
          key: 'type-n-enteroendocrine-cell',
        },
        {
          exp: '(?i)^enteroendocrine cell$',
          key: 'enteroendocrine-cell',
        },
        {
          exp: '(?i)^endocrine cell$',
          key: 'endocrine-cell',
        },
        {
          exp: '(?i)^unipolar brush cell$',
          key: 'unipolar-brush-cell',
        },
        {
          exp: '(?i)^brush cell$',
          key: 'brush-cell',
        },
        {
          exp: '(?i)^unswitched memory B cell$',
          key: 'unswitched-memory-b-cell',
        },
        {
          exp: '(?i)^memory B cell$',
          key: 'memory-b-cell',
        },
        {
          exp: '(?i)^B cell$',
          key: 'b-cell',
        },
        {
          exp: '(?i)^urethra urothelial cell$',
          key: 'urethra-urothelial-cell',
        },
        {
          exp: '(?i)^urothelial cell$',
          key: 'urothelial-cell',
        },
        {
          exp: '(?i)^uterine smooth muscle cell$',
          key: 'uterine-smooth-muscle-cell',
        },
        {
          exp: '(?i)^vasa recta ascending limb cell$',
          key: 'vasa-recta-ascending-limb-cell',
        },
        {
          exp: '(?i)^vasa recta descending limb cell$',
          key: 'vasa-recta-descending-limb-cell',
        },
        {
          exp: '(?i)^vascular associated smooth muscle cell$',
          key: 'vascular-associated-smooth-muscle-cell',
        },
        {
          exp: '(?i)^smooth muscle cell$',
          key: 'smooth-muscle-cell',
        },
        {
          exp: '(?i)^vascular leptomeningeal cell \\\\(Mmus\\\\)$',
          key: 'vascular-leptomeningeal-cell-mmus',
        },
        {
          exp: '(?i)^vascular leptomeningeal cell$',
          key: 'vascular-leptomeningeal-cell',
        },
        {
          exp: '(?i)^leptomeningeal cell$',
          key: 'leptomeningeal-cell',
        },
        {
          exp: '(?i)^vascular lymphangioblast$',
          key: 'vascular-lymphangioblast',
        },
        {
          exp: '(?i)^vein endothelial cell$',
          key: 'vein-endothelial-cell',
        },
        {
          exp: '(?i)^endothelial cell$',
          key: 'endothelial-cell',
        },
        {
          exp: '(?i)^ventricular cardiac muscle cell$',
          key: 'ventricular-cardiac-muscle-cell',
        },
        {
          exp: '(?i)^cardiac muscle cell$',
          key: 'cardiac-muscle-cell',
        },
        {
          exp: '(?i)^muscle cell$',
          key: 'muscle-cell',
        },
        {
          exp: '(?i)^cell$',
          key: 'cell',
        },
        {
          exp: '(?i)^vip GABAergic cortical interneuron$',
          key: 'vip-gabaergic-cortical-interneuron',
        },
        {
          exp: '(?i)^cortical interneuron$',
          key: 'cortical-interneuron',
        },
        {
          exp: '(?i)^interneuron$',
          key: 'interneuron',
        },
        {
          exp: '(?i)^visceromotor neuron$',
          key: 'visceromotor-neuron',
        },
        {
          exp: '(?i)^motor neuron$',
          key: 'motor-neuron',
        },
        {
          exp: '(?i)^neuron$',
          key: 'neuron',
        },
      ],
    },
    azimuth_coarse_cell_types: {
      match: [
        {
          exp: '^%s$' % std.strReplace(k, '+', '\\+'),
          key: std.strReplace(std.strReplace(std.strReplace(
            std.strReplace(std.asciiLower(k), '+ ', '-'), ' - ', '-'
          ), ' ', '-'), '/', '-'),
        }
        for k in [
          'Low UMI Barcode',
          'Coarse Not Confidently Annotated',
          'AT1',
          'AT2',
          'Adipocyte',
          'Amacrine cell',
          'Astrocyte',
          'Atrial cardiomyocyte',
          'BNC2+ neuron',
          'CLP',
          'Cone cell',
          'Cone-OFF-bipolar cell',
          'Cone-ON-bipolar cell',
          'Cycling T/NK cell',
          'Cycling epithelial cell',
          'Cycling immune cell',
          'Cycling myeloid cell',
          'Cycling perivascular cell',
          'Dendritic cell',
          'Doublet like cell',
          'ETV1+ neuron',
          'Early SOX4+ erythroblast',
          'Embryonic cell',
          'Endocrine cell',
          'Endothelial cell',
          'Epithelial cell of adrenal gland',
          'Epithelial cell of bladder',
          'Epithelial cell of breast',
          'Epithelial cell of eye',
          'Epithelial cell of gingiva',
          'Epithelial cell of intestine',
          'Epithelial cell of liver',
          'Epithelial cell of lung',
          'Epithelial cell of ovary',
          'Epithelial cell of pancreas',
          'Epithelial cell of skin',
          'Epithelial cell of stomach',
          'Epithelial cell of testis',
          'Epithelial cell of thymus',
          'Epithelial cell of tonsil',
          'Excitatory neuron',
          'Fibroblast',
          'Ganglion cell',
          'Germ cell',
          'Glial cell',
          'Granulocyte',
          'HSC',
          'Hematopoietic precursor cell',
          'Horizontal cell',
          'ILC',
          'Inhibitory neuron',
          'Intermediate EPCAM+ erythroblast',
          'KRT14 basal cell',
          'KRT14/ECM1 basal cell',
          'ECM1 basal cell',
          'Late hemoglobin+ erythroblast',
          'MEP',
          'Macrophage',
          'Mesothelial cell',
          'Miscellaneous neuron',
          'Monocyte/DC',
          'Monocyte',
          'NK cell',
          'Neuroblast',
          'Oligodendrocyte',
          'Osteoclast',
          'PCP4L1+ bipolar cell',
          'PDAC related epithelial cell',
          'Pericyte',
          'Platelet',
          'Podocyte',
          'PreB cell',
          'ProB cell',
          'B cell',
          'Red blood cell',
          'Renal epithelial cell - Loop of Henle',
          'Renal epithelial cell - distal tubules',
          'Renal epithelial cell - proximal tubules',
          'Rod cell',
          'Rod-bipolar cell',
          'Skeletal muscle cell',
          'T cell',
          'Type A intercalated cell',
          'Type B intercalated cell',
          'Unassigned',
          'Vascular smooth muscle cell',
          'Ventricular cardiomyocyte',
        ]
      ],
    },
    azimuth_fine_cell_types: {
      match: [
        {
          exp: '^Low UMI Barcode$',
          key: 'low-umi-barcode',
        },
        {
          exp: '^Fine Not Confidently Annotated$',
          key: 'fine-not-confidently-annotated',
        },
        {
          exp: '^ACTA2\\+ EC$',
          key: 'acta2-ec',
        },
        {
          exp: '^AFP\\+ epithelial cell$',
          key: 'afp-epithelial-cell',
        },
        {
          exp: '^ANGPTL7\\+NRP2\\+ fibroblast$',
          key: 'angptl7nrp2-fibroblast',
        },
        {
          exp: '^ASDC$',
          key: 'asdc',
        },
        {
          exp: '^AT1$',
          key: 'at1',
        },
        {
          exp: '^AT2$',
          key: 'at2',
        },
        {
          exp: '^Acinar cell$',
          key: 'acinar-cell',
        },
        {
          exp: '^Adipocyte$',
          key: 'adipocyte',
        },
        {
          exp: '^Adipose ITLN1 fibroblast$',
          key: 'adipose-itln1-fibroblast',
        },
        {
          exp: '^Aerocyte$',
          key: 'aerocyte',
        },
        {
          exp: '^Agonist T cell$',
          key: 'agonist-t-cell',
        },
        {
          exp: '^Airway basal cell$',
          key: 'airway-basal-cell',
        },
        {
          exp: '^Alpha cell$',
          key: 'alpha-cell',
        },
        {
          exp: '^Alveolar macrophage$',
          key: 'alveolar-macrophage',
        },
        {
          exp: '^Amacrine cell$',
          key: 'amacrine-cell',
        },
        {
          exp: '^Amygdala excitatory neuron$',
          key: 'amygdala-excitatory-neuron',
        },
        {
          exp: '^Arterial EC$',
          key: 'arterial-ec',
        },
        {
          exp: '^Astrocyte$',
          key: 'astrocyte',
        },
        {
          exp: '^Atrial cardiomyocyte$',
          key: 'atrial-cardiomyocyte',
        },
        {
          exp: '^BGN\\+ EC$',
          key: 'bgn-ec',
        },
        {
          exp: '^BNC2\\+ neuron$',
          key: 'bnc2-neuron',
        },
        {
          exp: '^BNC2\\+ZFPM2\\+ fibroblast$',
          key: 'bnc2zfpm2-fibroblast',
        },
        {
          exp: '^Basal cell of bladder$',
          key: 'basal-cell-of-bladder',
        },
        {
          exp: '^Basophil$',
          key: 'basophil',
        },
        {
          exp: '^Beta cell$',
          key: 'beta-cell',
        },
        {
          exp: '^C2orf40\\+ fibroblast$',
          key: 'c2orf40-fibroblast',
        },
        {
          exp: '^C3\\+ fibroblast$',
          key: 'c3-fibroblast',
        },
        {
          exp: '^CCL11\\+ fibroblast$',
          key: 'ccl11-fibroblast',
        },
        {
          exp: '^CCL2\\+SFRP2\\+ fibroblast$',
          key: 'ccl2sfrp2-fibroblast',
        },
        {
          exp: '^CD14 monocyte$',
          key: 'cd14-monocyte',
        },
        {
          exp: '^CD14/CD16 monocyte$',
          key: 'cd14-cd16-monocyte',
        },
        {
          exp: '^CD16 NK cell$',
          key: 'cd16-nk-cell',
        },
        {
          exp: '^CD16 monocyte$',
          key: 'cd16-monocyte',
        },
        {
          exp: '^CD55\\+SEMA3C\\+ fibroblast$',
          key: 'cd55sema3c-fibroblast',
        },
        {
          exp: '^CD56 NK cell$',
          key: 'cd56-nk-cell',
        },
        {
          exp: '^CD8aa\\+ T cell$',
          key: 'cd8aa-t-cell',
        },
        {
          exp: '^CDCA7\\+ fibroblast$',
          key: 'cdca7-fibroblast',
        },
        {
          exp: '^CFD\\+MGP\\+ fibroblast$',
          key: 'cfdmgp-fibroblast',
        },
        {
          exp: '^CLDN3\\+ epithelial cell$',
          key: 'cldn3-epithelial-cell',
        },
        {
          exp: '^CLP$',
          key: 'clp',
        },
        {
          exp: '^CLU\\+ fibroblast$',
          key: 'clu-fibroblast',
        },
        {
          exp: '^CRABP1\\+ mesothelial cell$',
          key: 'crabp1-mesothelial-cell',
        },
        {
          exp: '^CXCL\\+ fibroblast$',
          key: 'cxcl-fibroblast',
        },
        {
          exp: '^CXCL13 exhausted CD8 T cell$',
          key: 'cxcl13-exhausted-cd8-t-cell',
        },
        {
          exp: '^CXCL14 inhibitory neuron$',
          key: 'cxcl14-inhibitory-neuron',
        },
        {
          exp: '^Capillary EC$',
          key: 'capillary-ec',
        },
        {
          exp: '^Capsular cell$',
          key: 'capsular-cell',
        },
        {
          exp: '^Chief cell$',
          key: 'chief-cell',
        },
        {
          exp: '^Chromaffin cell$',
          key: 'chromaffin-cell',
        },
        {
          exp: '^Ciliary body cell$',
          key: 'ciliary-body-cell',
        },
        {
          exp: '^Club cell$',
          key: 'club-cell',
        },
        {
          exp: '^Colonocyte$',
          key: 'colonocyte',
        },
        {
          exp: '^Cone cell$',
          key: 'cone-cell',
        },
        {
          exp: '^Cone-OFF-bipolar cell$',
          key: 'cone-off-bipolar-cell',
        },
        {
          exp: '^Cone-ON-bipolar cell$',
          key: 'cone-on-bipolar-cell',
        },
        {
          exp: '^Corneal epithelial cell$',
          key: 'corneal-epithelial-cell',
        },
        {
          exp: '^Cortical thymic epithelial cell$',
          key: 'cortical-thymic-epithelial-cell',
        },
        {
          exp: '^Cycling B$',
          key: 'cycling-b',
        },
        {
          exp: '^Cycling DN/DP T cell$',
          key: 'cycling-dn-dp-t-cell',
        },
        {
          exp: '^Cycling EC$',
          key: 'cycling-ec',
        },
        {
          exp: '^Cycling PDAC specific ductal cell$',
          key: 'cycling-pdac-specific-ductal-cell',
        },
        {
          exp: '^Cycling T cell$',
          key: 'cycling-t-cell',
        },
        {
          exp: '^Cycling T/NK cell$',
          key: 'cycling-t-nk-cell',
        },
        {
          exp: '^Cycling epithelial cell$',
          key: 'cycling-epithelial-cell',
        },
        {
          exp: '^Cycling fibroblast$',
          key: 'cycling-fibroblast',
        },
        {
          exp: '^Cycling immune cell$',
          key: 'cycling-immune-cell',
        },
        {
          exp: '^Cycling macrophage$',
          key: 'cycling-macrophage',
        },
        {
          exp: '^Cycling mast cell$',
          key: 'cycling-mast-cell',
        },
        {
          exp: '^Cycling mesodermal cell$',
          key: 'cycling-mesodermal-cell',
        },
        {
          exp: '^Cycling mesothelial cell$',
          key: 'cycling-mesothelial-cell',
        },
        {
          exp: '^Cycling monocyte$',
          key: 'cycling-monocyte',
        },
        {
          exp: '^Cycling perivascular cell$',
          key: 'cycling-perivascular-cell',
        },
        {
          exp: '^Cycling thymic epithelial cell$',
          key: 'cycling-thymic-epithelial-cell',
        },
        {
          exp: '^Cytotrophoblast$',
          key: 'cytotrophoblast',
        },
        {
          exp: '^D1 medium spiny neuron$',
          key: 'd1-medium-spiny-neuron',
        },
        {
          exp: '^D2 medium spiny neuron$',
          key: 'd2-medium-spiny-neuron',
        },
        {
          exp: '^DIAPH3\\+ macrophage$',
          key: 'diaph3-macrophage',
        },
        {
          exp: '^Delta cell$',
          key: 'delta-cell',
        },
        {
          exp: '^Double negative T cell$',
          key: 'double-negative-t-cell',
        },
        {
          exp: '^Double positive T cell$',
          key: 'double-positive-t-cell',
        },
        {
          exp: '^Ductal/EC doublet like cell$',
          key: 'ductal-ec-doublet-like-cell',
        },
        {
          exp: '^ETV1\\+ neuron$',
          key: 'etv1-neuron',
        },
        {
          exp: '^Early SOX4\\+ erythroblast$',
          key: 'early-sox4-erythroblast',
        },
        {
          exp: '^Eccentric medium spiny neuron$',
          key: 'eccentric-medium-spiny-neuron',
        },
        {
          exp: '^Endocardial cell$',
          key: 'endocardial-cell',
        },
        {
          exp: '^Endocrine cell$',
          key: 'endocrine-cell',
        },
        {
          exp: '^Endothelial‑mesenchymal transition EC$',
          key: 'endothelial-mesenchymal-transition-ec',
        },
        {
          exp: '^Enterochromaffin cell$',
          key: 'enterochromaffin-cell',
        },
        {
          exp: '^Enterocyte$',
          key: 'enterocyte',
        },
        {
          exp: '^Enteroendocrine cell$',
          key: 'enteroendocrine-cell',
        },
        {
          exp: '^EoBaMa Progenitor$',
          key: 'eobama-progenitor',
        },
        {
          exp: '^Eosinophil$',
          key: 'eosinophil',
        },
        {
          exp: '^Epithelial cell of conjunctiva$',
          key: 'epithelial-cell-of-conjunctiva',
        },
        {
          exp: '^Epithelial cell of tonsil$',
          key: 'epithelial-cell-of-tonsil',
        },
        {
          exp: '^Epsilon cell$',
          key: 'epsilon-cell',
        },
        {
          exp: '^Extravillous trophoblast$',
          key: 'extravillous-trophoblast',
        },
        {
          exp: '^Eye melanocyte$',
          key: 'eye-melanocyte',
        },
        {
          exp: '^F3\\+ fibroblast$',
          key: 'f3-fibroblast',
        },
        {
          exp: '^FBLN2\\+APOD\\+ fibroblast$',
          key: 'fbln2apod-fibroblast',
        },
        {
          exp: '^FRZB\\+ fibroblast$',
          key: 'frzb-fibroblast',
        },
        {
          exp: '^Fetal HBG\\+ erythrocyte$',
          key: 'fetal-hbg-erythrocyte',
        },
        {
          exp: '^Fetal KCNQ1OT1\\+MIR503HG\\+ fibroblast$',
          key: 'fetal-kcnq1ot1mir503hg-fibroblast',
        },
        {
          exp: '^Fibroblast/Granulosa doublet like cell$',
          key: 'fibroblast-granulosa-doublet-like-cell',
        },
        {
          exp: '^Follicular DC$',
          key: 'follicular-dc',
        },
        {
          exp: '^G0S2\\+PPP1R14A\\+ fibroblast$',
          key: 'g0s2ppp1r14a-fibroblast',
        },
        {
          exp: '^G2/M phase myeloid cell$',
          key: 'g2-m-phase-myeloid-cell',
        },
        {
          exp: '^GABRP\\+ epithelial cell$',
          key: 'gabrp-epithelial-cell',
        },
        {
          exp: '^GATA4\\+ mesothelial cell$',
          key: 'gata4-mesothelial-cell',
        },
        {
          exp: '^GMP$',
          key: 'gmp',
        },
        {
          exp: '^GREB1\\+ fibroblast$',
          key: 'greb1-fibroblast',
        },
        {
          exp: '^GZMB CD8 T cell$',
          key: 'gzmb-cd8-t-cell',
        },
        {
          exp: '^GZMK CD8 T cell$',
          key: 'gzmk-cd8-t-cell',
        },
        {
          exp: '^GZMK cytotoxic CD4 T cell$',
          key: 'gzmk-cytotoxic-cd4-t-cell',
        },
        {
          exp: '^Ganglion cell$',
          key: 'ganglion-cell',
        },
        {
          exp: '^Germinal center B cell$',
          key: 'germinal-center-b-cell',
        },
        {
          exp: '^Glial cell of intestine$',
          key: 'glial-cell-of-intestine',
        },
        {
          exp: '^Goblet cell$',
          key: 'goblet-cell',
        },
        {
          exp: '^Granulosa cell$',
          key: 'granulosa-cell',
        },
        {
          exp: '^HSC$',
          key: 'hsc',
        },
        {
          exp: '^Hepatocyte$',
          key: 'hepatocyte',
        },
        {
          exp: '^Hippocampal CA1-3 neuron$',
          key: 'hippocampal-ca1-3-neuron',
        },
        {
          exp: '^Hippocampal CA4 neuron$',
          key: 'hippocampal-ca4-neuron',
        },
        {
          exp: '^Hippocampal dentate gyrus neuron$',
          key: 'hippocampal-dentate-gyrus-neuron',
        },
        {
          exp: '^Hofbauer cell$',
          key: 'hofbauer-cell',
        },
        {
          exp: '^Horizontal cell$',
          key: 'horizontal-cell',
        },
        {
          exp: '^IGFBP1\\+ fibroblast$',
          key: 'igfbp1-fibroblast',
        },
        {
          exp: '^IGFBP6\\+APOD\\+ fibroblast$',
          key: 'igfbp6apod-fibroblast',
        },
        {
          exp: '^IGFBP6\\+SFRP4\\+ fibroblast$',
          key: 'igfbp6sfrp4-fibroblast',
        },
        {
          exp: '^IL1B\\+ macrophage$',
          key: 'il1b-macrophage',
        },
        {
          exp: '^ILC$',
          key: 'ilc',
        },
        {
          exp: '^INF-activated CD8 T cell$',
          key: 'inf-activated-cd8-t-cell',
        },
        {
          exp: '^INF-activated T cell$',
          key: 'inf-activated-t-cell',
        },
        {
          exp: '^INF-activated macrophage$',
          key: 'inf-activated-macrophage',
        },
        {
          exp: '^INF-activated trophoblast$',
          key: 'inf-activated-trophoblast',
        },
        {
          exp: '^Intermediate EPCAM\\+ erythroblast$',
          key: 'intermediate-epcam-erythroblast',
        },
        {
          exp: '^Intestinal stem cell$',
          key: 'intestinal-stem-cell',
        },
        {
          exp: '^KCNN3\\+ fibroblast$',
          key: 'kcnn3-fibroblast',
        },
        {
          exp: '^KLRB1 CD8 T cell$',
          key: 'klrb1-cd8-t-cell',
        },
        {
          exp: '^KLRB1 cytotoxic CD4 T$',
          key: 'klrb1-cytotoxic-cd4-t',
        },
        {
          exp: '^KRT14 basal cell$',
          key: 'krt14-basal-cell',
        },
        {
          exp: '^KRT14/ECM1 basal cell$',
          key: 'krt14-ecm1-basal-cell',
        },
        {
          exp: '^ECM1 basal cell$',
          key: 'ecm1-basal-cell',
        },
        {
          exp: '^Keratinocyte$',
          key: 'keratinocyte',
        },
        {
          exp: '^L2/3 excitatory neuron$',
          key: 'l2-3-excitatory-neuron',
        },
        {
          exp: '^L4 excitatory neuron$',
          key: 'l4-excitatory-neuron',
        },
        {
          exp: '^L5 excitatory neuron$',
          key: 'l5-excitatory-neuron',
        },
        {
          exp: '^L5/6 excitatory neuron$',
          key: 'l5-6-excitatory-neuron',
        },
        {
          exp: '^L5b excitatory neuron$',
          key: 'l5b-excitatory-neuron',
        },
        {
          exp: '^L6 excitatory neuron$',
          key: 'l6-excitatory-neuron',
        },
        {
          exp: '^LMPP$',
          key: 'lmpp',
        },
        {
          exp: '^LYVE1 macrophage$',
          key: 'lyve1-macrophage',
        },
        {
          exp: '^Lactocyte$',
          key: 'lactocyte',
        },
        {
          exp: '^Late hemoglobin\\+ erythroblast$',
          key: 'late-hemoglobin-erythroblast',
        },
        {
          exp: '^Leydig cell$',
          key: 'leydig-cell',
        },
        {
          exp: '^Lower rhombic lip-derived neuron$',
          key: 'lower-rhombic-lip-derived-neuron',
        },
        {
          exp: '^Luminal cell$',
          key: 'luminal-cell',
        },
        {
          exp: '^Lymphatic EC$',
          key: 'lymphatic-ec',
        },
        {
          exp: '^M1 macrophage$',
          key: 'm1-macrophage',
        },
        {
          exp: '^MAIT cell$',
          key: 'mait-cell',
        },
        {
          exp: '^MDK\\+NREP\\+ fibroblast$',
          key: 'mdknrep-fibroblast',
        },
        {
          exp: '^MEP$',
          key: 'mep',
        },
        {
          exp: '^MMP1\\+ fibroblast$',
          key: 'mmp1-fibroblast',
        },
        {
          exp: '^MMP19 macrophage$',
          key: 'mmp19-macrophage',
        },
        {
          exp: '^MPP$',
          key: 'mpp',
        },
        {
          exp: '^MTRNR2L1 B cell$',
          key: 'mtrnr2l1-b-cell',
        },
        {
          exp: '^Mammary basal cell$',
          key: 'mammary-basal-cell',
        },
        {
          exp: '^Mammary luminal cell$',
          key: 'mammary-luminal-cell',
        },
        {
          exp: '^Mammary luminal progenitor$',
          key: 'mammary-luminal-progenitor',
        },
        {
          exp: '^Mammillary body neuron$',
          key: 'mammillary-body-neuron',
        },
        {
          exp: '^Mast cell$',
          key: 'mast-cell',
        },
        {
          exp: '^Medullary thymic epithelial cell$',
          key: 'medullary-thymic-epithelial-cell',
        },
        {
          exp: '^Melanocyte$',
          key: 'melanocyte',
        },
        {
          exp: '^Memory B cell precursor$',
          key: 'memory-b-cell-precursor',
        },
        {
          exp: '^Memory B cell$',
          key: 'memory-b-cell',
        },
        {
          exp: '^Memory CD4 T cell$',
          key: 'memory-cd4-t-cell',
        },
        {
          exp: '^Memory T cell$',
          key: 'memory-t-cell',
        },
        {
          exp: '^Merkel cell carcinoma enriched epithelial cell$',
          key: 'merkel-cell-carcinoma-enriched-epithelial-cell',
        },
        {
          exp: '^Merkel cell$',
          key: 'merkel-cell',
        },
        {
          exp: '^Mesenchymal stromal cell$',
          key: 'mesenchymal-stromal-cell',
        },
        {
          exp: '^Metallothionein\\+ epithelial cell$',
          key: 'metallothionein-epithelial-cell',
        },
        {
          exp: '^Metallothionein\\+ fibroblast$',
          key: 'metallothionein-fibroblast',
        },
        {
          exp: '^Metallothionein\\+ macrophage$',
          key: 'metallothionein-macrophage',
        },
        {
          exp: '^Microfold cell$',
          key: 'microfold-cell',
        },
        {
          exp: '^Microglia$',
          key: 'microglia',
        },
        {
          exp: '^Midbrain-derived inhibitory neuron$',
          key: 'midbrain-derived-inhibitory-neuron',
        },
        {
          exp: '^MkP$',
          key: 'mkp',
        },
        {
          exp: '^Monocyte/DC$',
          key: 'monocyte-dc',
        },
        {
          exp: '^Muller glia$',
          key: 'muller-glia',
        },
        {
          exp: '^Multiciliated cell$',
          key: 'multiciliated-cell',
        },
        {
          exp: '^Myofibroblast$',
          key: 'myofibroblast',
        },
        {
          exp: '^NR4A3 fibroblast$',
          key: 'nr4a3-fibroblast',
        },
        {
          exp: '^Naive B cell$',
          key: 'naive-b-cell',
        },
        {
          exp: '^Naive CD4 T cell$',
          key: 'naive-cd4-t-cell',
        },
        {
          exp: '^Naive CD8 T cell$',
          key: 'naive-cd8-t-cell',
        },
        {
          exp: '^Naive T cell$',
          key: 'naive-t-cell',
        },
        {
          exp: '^Neuroblast$',
          key: 'neuroblast',
        },
        {
          exp: '^Neuroendocrine cell$',
          key: 'neuroendocrine-cell',
        },
        {
          exp: '^Neutrophil precursor$',
          key: 'neutrophil-precursor',
        },
        {
          exp: '^Neutrophil$',
          key: 'neutrophil',
        },
        {
          exp: '^Oligodendrocyte$',
          key: 'oligodendrocyte',
        },
        {
          exp: '^Oocyte$',
          key: 'oocyte',
        },
        {
          exp: '^Osteoclast$',
          key: 'osteoclast',
        },
        {
          exp: '^Ovarian cancer related epithelial cell$',
          key: 'ovarian-cancer-related-epithelial-cell',
        },
        {
          exp: '^PAX2 inhibitory neuron$',
          key: 'pax2-inhibitory-neuron',
        },
        {
          exp: '^PCOLCE2 fibroblast$',
          key: 'pcolce2-fibroblast',
        },
        {
          exp: '^PCP4L1\\+ bipolar cell$',
          key: 'pcp4l1-bipolar-cell',
        },
        {
          exp: '^PDAC specific ductal cell$',
          key: 'pdac-specific-ductal-cell',
        },
        {
          exp: '^PLA2G2A\\+IGFBP4\\+APOC1\\+ fibroblast$',
          key: 'pla2g2aigfbp4apoc1-fibroblast',
        },
        {
          exp: '^PLAC9\\+ fibroblast$',
          key: 'plac9-fibroblast',
        },
        {
          exp: '^PLCH1 L4/5 excitatory neuron$',
          key: 'plch1-l4-5-excitatory-neuron',
        },
        {
          exp: '^POSTN fibroblast$',
          key: 'postn-fibroblast',
        },
        {
          exp: '^PP cell$',
          key: 'pp-cell',
        },
        {
          exp: '^PVALB inhibitory neuron$',
          key: 'pvalb-inhibitory-neuron',
        },
        {
          exp: '^Pancreatic ductal cell$',
          key: 'pancreatic-ductal-cell',
        },
        {
          exp: '^Paneth cell$',
          key: 'paneth-cell',
        },
        {
          exp: '^Pericyte$',
          key: 'pericyte',
        },
        {
          exp: '^Peritubular myoid cell$',
          key: 'peritubular-myoid-cell',
        },
        {
          exp: '^Placenta CAMP\\+ monocyte$',
          key: 'placenta-camp-monocyte',
        },
        {
          exp: '^Placenta defensin\\+ monocyte$',
          key: 'placenta-defensin-monocyte',
        },
        {
          exp: '^Plasma cell$',
          key: 'plasma-cell',
        },
        {
          exp: '^Platelet$',
          key: 'platelet',
        },
        {
          exp: '^Podocyte$',
          key: 'podocyte',
        },
        {
          exp: '^Pre-ProB cell$',
          key: 'pre-prob-cell',
        },
        {
          exp: '^Pre-cDC$',
          key: 'pre-cdc',
        },
        {
          exp: '^Pre-pDC$',
          key: 'pre-pdc',
        },
        {
          exp: '^PreB cell$',
          key: 'preb-cell',
        },
        {
          exp: '^Primordial germ cell$',
          key: 'primordial-germ-cell',
        },
        {
          exp: '^ProB cell$',
          key: 'prob-cell',
        },
        {
          exp: '^Pulmonary ionocyte$',
          key: 'pulmonary-ionocyte',
        },
        {
          exp: '^Pyramidal neuron$',
          key: 'pyramidal-neuron',
        },
        {
          exp: '^RNASE1\\+ fibroblast$',
          key: 'rnase1-fibroblast',
        },
        {
          exp: '^Renal epithelial cell - Loop of Henle$',
          key: 'renal-epithelial-cell-loop-of-henle',
        },
        {
          exp: '^Renal epithelial cell - distal tubules$',
          key: 'renal-epithelial-cell-distal-tubules',
        },
        {
          exp: '^Renal epithelial cell - proximal tubules$',
          key: 'renal-epithelial-cell-proximal-tubules',
        },
        {
          exp: '^Retinal pigment epithelial cell$',
          key: 'retinal-pigment-epithelial-cell',
        },
        {
          exp: '^Rod cell$',
          key: 'rod-cell',
        },
        {
          exp: '^Rod-bipolar cell$',
          key: 'rod-bipolar-cell',
        },
        {
          exp: '^S phase myeloid cell$',
          key: 's-phase-myeloid-cell',
        },
        {
          exp: '^SCGB2A1\\+ endocrine cell$',
          key: 'scgb2a1-endocrine-cell',
        },
        {
          exp: '^SCGN inhibitory neuron$',
          key: 'scgn-inhibitory-neuron',
        },
        {
          exp: '^SCN7A fibroblast$',
          key: 'scn7a-fibroblast',
        },
        {
          exp: '^SPARCL1\\+ fibroblast$',
          key: 'sparcl1-fibroblast',
        },
        {
          exp: '^SPP1 macrophage$',
          key: 'spp1-macrophage',
        },
        {
          exp: '^SST inhibitory neuron$',
          key: 'sst-inhibitory-neuron',
        },
        {
          exp: '^SV2C inhibitory neuron$',
          key: 'sv2c-inhibitory-neuron',
        },
        {
          exp: '^Schwann cell$',
          key: 'schwann-cell',
        },
        {
          exp: '^Sertoli and Interstitial cell common progenitor$',
          key: 'sertoli-and-interstitial-cell-common-progenitor',
        },
        {
          exp: '^Sertoli cell$',
          key: 'sertoli-cell',
        },
        {
          exp: '^Skeletal muscle cell$',
          key: 'skeletal-muscle-cell',
        },
        {
          exp: '^Spermatocyte$',
          key: 'spermatocyte',
        },
        {
          exp: '^Splatter neuron$',
          key: 'splatter-neuron',
        },
        {
          exp: '^Squamous epithelial cell$',
          key: 'squamous-epithelial-cell',
        },
        {
          exp: '^Steroidogenic cell$',
          key: 'steroidogenic-cell',
        },
        {
          exp: '^Supporting-like PAX8\\+ cell$',
          key: 'supporting-like-pax8-cell',
        },
        {
          exp: '^Surface mucous cell$',
          key: 'surface-mucous-cell',
        },
        {
          exp: '^Sympathoadrenal cell$',
          key: 'sympathoadrenal-cell',
        },
        {
          exp: '^Sympathoblast$',
          key: 'sympathoblast',
        },
        {
          exp: '^Syncytiotrophoblast$',
          key: 'syncytiotrophoblast',
        },
        {
          exp: '^TCIM\\+NRG1\\+ fibroblast$',
          key: 'tcimnrg1-fibroblast',
        },
        {
          exp: '^THY1\\+ fibroblast$',
          key: 'thy1-fibroblast',
        },
        {
          exp: '^TNC fibroblast$',
          key: 'tnc-fibroblast',
        },
        {
          exp: '^TSHZ2 L4/5 excitatory neuron$',
          key: 'tshz2-l4-5-excitatory-neuron',
        },
        {
          exp: '^Tfh cell$',
          key: 'tfh-cell',
        },
        {
          exp: '^Thalamic excitatory neuron$',
          key: 'thalamic-excitatory-neuron',
        },
        {
          exp: '^Thymic ciliated cell$',
          key: 'thymic-ciliated-cell',
        },
        {
          exp: '^Thymic myoid cell$',
          key: 'thymic-myoid-cell',
        },
        {
          exp: '^Tip EC$',
          key: 'tip-ec',
        },
        {
          exp: '^Tissue-resident NK cell$',
          key: 'tissue-resident-nk-cell',
        },
        {
          exp: '^Tonsillar crypt epithelial cell$',
          key: 'tonsillar-crypt-epithelial-cell',
        },
        {
          exp: '^Transit-amplifying cell$',
          key: 'transit-amplifying-cell',
        },
        {
          exp: '^Treg cell$',
          key: 'treg-cell',
        },
        {
          exp: '^Tuft cell$',
          key: 'tuft-cell',
        },
        {
          exp: '^Type A intercalated cell$',
          key: 'type-a-intercalated-cell',
        },
        {
          exp: '^Type B intercalated cell$',
          key: 'type-b-intercalated-cell',
        },
        {
          exp: '^Unassigned$',
          key: 'unassigned',
        },
        {
          exp: '^Upper rhombic lip-derived neuron$',
          key: 'upper-rhombic-lip-derived-neuron',
        },
        {
          exp: '^Uveal melanocyte$',
          key: 'uveal-melanocyte',
        },
        {
          exp: '^VIP inhibitory neuron$',
          key: 'vip-inhibitory-neuron',
        },
        {
          exp: '^Vascular smooth muscle cell$',
          key: 'vascular-smooth-muscle-cell',
        },
        {
          exp: '^Venous EC$',
          key: 'venous-ec',
        },
        {
          exp: '^Ventricular cardiomyocyte$',
          key: 'ventricular-cardiomyocyte',
        },
        {
          exp: '^WISP2\\+ fibroblast$',
          key: 'wisp2-fibroblast',
        },
        {
          exp: '^XCL1 NK cell$',
          key: 'xcl1-nk-cell',
        },
        {
          exp: '^ZEB2\\+ mesoderm$',
          key: 'zeb2-mesoderm',
        },
        {
          exp: '^abT \\(entry\\) cell$',
          key: 'abt-entry-cell',
        },
        {
          exp: '^cDC1$',
          key: 'cdc1',
        },
        {
          exp: '^cDC2$',
          key: 'cdc2',
        },
        {
          exp: '^cholangiocyte$',
          key: 'cholangiocyte',
        },
        {
          exp: '^gdT cell$',
          key: 'gdt-cell',
        },
        {
          exp: '^mregDC$',
          key: 'mregdc',
        },
        {
          exp: '^pDC$',
          key: 'pdc',
        },
      ],
    },
  },

  // Presets shared by Cell Ranger and Space Ranger, but not ARC.
  cr_sr_presets: {
    reference_genomes: {
      match: [
        {
          exp: '(?i)(human)|(grch38)|(grch37)|(hu.custom)|(homo sapiens)|(hg38)|(homo_sapiens)|(homo.sapiens)|(hs38)|(hsapiens)|(h38)$',
          key: 'Human',
        },
        {
          exp: '(?i)(mouse)|(mm10)|(mm39)|(grcm39)|(mmus)|(musculus)|(mus)|(grcm)|(mma2024)|(mm_egfp)|(mice)|(mma2020)|(m39)|(mm38)|(pman)|(spret_eij)|(mm_sp)|(m10)$',
          key: 'Mouse',
        },
        {
          exp: '(?i)(rat)|(mrat)|(rn)|(rn6)|(rno6)$',
          key: 'Rat',
        },
        {
          exp: '(?i)(mmul)|(cyno)|(jacchus)|(macaca)|(mcalja1)|(mfas)|(monkey)|(macaque)|(rhesus)|(callithrix)|(marmoset)|(rhemac)$',
          key: 'Macaque',
        },
        {
          exp: '(?i)(danio)|(drerio)|(grcz11)|(danrer)|(zebrafish)|(d_rerio)$',
          key: 'Zebrafish',
        },
        {
          exp: '(?i)(agambiae)$',
          key: 'Mosquito',
        },
        {
          exp: '(?i)(pf3d7)|(hs_mng_pf)|(plasmodb)|(malaria)|(pfgab4)|(plasmodium)|(pfgab4)|(pfgb4)$',
          key: 'Malaria',
        },
        {
          exp: '(?i)(rabbit)|(ycun)|(yctolagus.cuniculus)$',
          key: 'Rabbit',
        },
        {
          exp: '(?i)(setaria)|(arabidopsis)|(caijin_osativa)|(arahy)|(tair10)|(zeamays)|(cowpea_v1.1)|(medicago_truncatula)$',
          key: 'Plant',
        },
        {
          exp: '(?i)(canine)|(cfam)|(canfam)|(canis.lupus)|(canis_lupus)|(dog)$',
          key: 'Dog',
        },
        {
          exp: '(?i)(sus_scrofa)|(swine)|(sscrofa)|(sscofa)|(pig)$',
          key: 'Pig',
        },
        {
          exp: '(?i)(dmel)|(drosophila)|(melanogaster)|(dm6)$',
          key: 'Fly',
        },
        {
          exp: '(?i)(elegans)$',
          key: 'Elegans',
        },
        {
          exp: '(?i)(amexg)|(ddv6)|(amexicanus)|(dd_chrfl)$',
          key: 'Axolotl',
        },
        {
          exp: '(?i)(gallus)|(grcg7b)|(leghn)|(galgal)|(chicken)$',
          key: 'Chicken',
        },
        {
          exp: '(?i)(ovis_aries)|(sheep)|(oarv2)|(ars-ui_ramb)$',
          key: 'Sheep',
        },
        {
          exp: '(?i)(clam)$',
          key: 'Clam',
        },
        {
          exp: '(?i)(aeruginosa)|(bacteria)$',
          key: 'Bacteria',
        },
        {
          exp: '(?i)(bos_taurus)|(btaurus)|(bosTau)|(bostau)|(bos-taurus)|(cow_)|(ars-ucd1)$',
          key: 'Cow',
        },
        {
          exp: '(?i)(xenopus_tropicalis)|(xenopus)$',
          key: 'Frog',
        },
        {
          exp: '(?i)(c_albicans)|(yeast)|(cerevisiae)|(kmar1.0)$',
          key: 'Yeast',
        },
        {
          exp: '(?i)(catfish)$',
          key: 'Catfish',
        },
        {
          exp: '(?i)(hse)|(equcab)|(horse)$',
          key: 'Horse',
        },
        {
          exp: '(?i)(crovir3)|(oyster)$',
          key: 'Oyster',
        },
        {
          exp: '(?i)(gondii)|(tgme49)|(toxoplasma)$',
          key: 'Toxoplasma',
        },
        {
          exp: '(?i)(mpantro3)|(chimpanzee)$',
          key: 'Chimpanzee',
        },
        {
          exp: '(?i)(rhifer)|(artibeus_jamaicensis)$',
          key: 'Bat',
        },
        {
          exp: '(?i)(hxb2)|(hiv)|(sivmac)$',
          key: 'HIV',
        },
        {
          exp: '(?i)(fish)$',
          key: 'Fish',
        },
        {
          exp: '(?i)(cviridis)|(snake)$',
          key: 'Snake',
        },
        {
          exp: '(?i)(kpetmar)|(lamprey)$',
          key: 'Lamprey',
        },
        {
          exp: '(?i)(maur)|(hamster)$',
          key: 'Hamster',
        },
        {
          exp: '(?i)(tursiops)|(dolphin)$',
          key: 'Dolphin',
        },
      ],
    },
    sequencer_platform: {
      match: [
        { exp: '^ilmn$', key: 'ilmn' },
        { exp: '^ont$', key: 'ont' },
        { exp: '^pb$', key: 'pb' },
        { exp: '^ult$', key: 'ult' },
        { exp: '^sbx$', key: 'sbx' },
        { exp: '.*', key: 'unknown' },
      ],
    },
    target_panel_hash: {
      match: [
        {
          key: 'Chromium Cynomolgus Monkey Transcriptome Probe Set v1.0.0',
          exp: '6de931d7dd5ff476b60ba54a9ee3275491c7ff15',
        },
        {
          key: 'Chromium Dog Transcriptome Probe Set v1.0.0',
          exp: 'd8876d2f5f2b34ec89f39c036f9a507d011c0601',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v1.0.1',
          exp: '3946c43ff94470bce5994d6c50ff9257ef122650',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v1.0.2',
          exp: '2e98f9dd1481353be2f31c8747c300f58e560d35',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set',
          exp: 'da653c08fca0c5c0e6cab30ddfe099c38d846a42',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v1.1.0',
          exp: 'cba2412e918a25f6f17ff1dede609e2bc58ca87c',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v2.0.0',
          exp: '7bdbe5cfbda2082f35b993f1ed5c59074576e0b2',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.0.1',
          exp: 'e54e103948ed4ff4fb868dcc116a7abb0193308d',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.0.2',
          exp: 'f8d3cf1a3ea9c982a78b5df9c3d94b9ec2817eb4',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set',
          exp: '7388d1718f8b62114816d6815202a6573b7aeca1',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.1.0',
          exp: '8c94cf42998dd9f8b21ad85b150dfa2e521e30e6',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.1.1',
          exp: '9be43d32826a017fd53069df576b83c0ab93bad7',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v2.0.0',
          exp: '93111cdbc73f9ed4d8a54802b2c77ab43c7b0825',
        },
        {
          key: 'Chromium Rat Transcriptome Probe Set v2.0.0',
          exp: 'b27fab796a7dbcbfcf0c8d9273b8828fa3c944cf',
        },
        {
          key: 'Chromium Rhesus Macaque Transcriptome Probe Set v1.0.0',
          exp: '0c640071fea2261113e83387bf898c8e5ae59cdb',
        },
        {
          key: 'Chromium and Visium Rat Transcriptome Probe Set v1.0.0',
          exp: '48b317b9e39be7b3b45c8cc5048ebb7ad019cc75',
        },
        {
          key: 'Visium Human Transcriptome Probe Set',
          exp: '3468b7772015471c99a8dc39af971d720b012869',
        },
        {
          key: 'Visium Human Transcriptome Probe Set v2.0',
          exp: 'dcd517e39ddd7e91c9e58caebf8451f69d62e38a',
        },
        {
          key: 'Visium Human Transcriptome Probe Set v2.1.0',
          exp: '900a8670482d7cdb02270f42ef21d267017bc200',
        },
        {
          key: 'Visium Mouse Transcriptome Probe Set',
          exp: '10426d8734e97c23fe1d5b72bcc1cb16885b5699',
        },
        {
          key: 'Visium Mouse Transcriptome Probe Set v2.0',
          exp: '45d58b6b029b088cf9e5b767c8941bd0a7aba297',
        },
        {
          key: 'Visium Mouse Transcriptome Probe Set v2.1.0',
          exp: '74a37b0d0a4af5f82473507f4760d0e0c78346ef',
        },
        {
          exp: '.*',
          key: 'unknown',
        },
      ],
    },
    target_panel_name: {
      match: [
        {
          key: 'Chromium Cynomolgus Monkey Transcriptome Probe Set v1.0.0',
          exp: '^Chromium Cynomolgus Monkey Transcriptome Probe Set v1.0.0$',
        },
        {
          key: 'Chromium Dog Transcriptome Probe Set v1.0.0',
          exp: '^Chromium Dog Transcriptome Probe Set v1.0.0$',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v1.0.1',
          exp: '^Chromium Human Transcriptome Probe Set v1.0.1$',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v1.0.2',
          exp: '^Chromium Human Transcriptome Probe Set v1.0.2$',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set',
          exp: '^Chromium Human Transcriptome Probe Set$',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v1.1.0',
          exp: '^Chromium Human Transcriptome Probe Set v1.1.0$',
        },
        {
          key: 'Chromium Human Transcriptome Probe Set v2.0.0',
          exp: '^Chromium Human Transcriptome Probe Set v2.0.0$',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.0.1',
          exp: '^Chromium Mouse Transcriptome Probe Set v1.0.1$',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.0.2',
          exp: '^Chromium Mouse Transcriptome Probe Set v1.0.2$',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set',
          exp: '^Chromium Mouse Transcriptome Probe Set$',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.1.0',
          exp: '^Chromium Mouse Transcriptome Probe Set v1.1.0$',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v1.1.1',
          exp: '^Chromium Mouse Transcriptome Probe Set v1.1.1$',
        },
        {
          key: 'Chromium Mouse Transcriptome Probe Set v2.0.0',
          exp: '^Chromium Mouse Transcriptome Probe Set v2.0.0$',
        },
        {
          key: 'Chromium Rat Transcriptome Probe Set v2.0.0',
          exp: '^Chromium Rat Transcriptome Probe Set v2.0.0$',
        },
        {
          key: 'Chromium Rhesus Macaque Transcriptome Probe Set v1.0.0',
          exp: '^Chromium Rhesus Macaque Transcriptome Probe Set v1.0.0$',
        },
        {
          key: 'Chromium and Visium Rat Transcriptome Probe Set v1.0.0',
          exp: '^Chromium and Visium Rat Transcriptome Probe Set v1.0.0$',
        },
        {
          key: 'Visium Human Transcriptome Probe Set',
          exp: '^Visium Human Transcriptome Probe Set$',
        },
        {
          key: 'Visium Human Transcriptome Probe Set v2.0',
          exp: '^Visium Human Transcriptome Probe Set v2.0$',
        },
        {
          key: 'Visium Human Transcriptome Probe Set v2.1.0',
          exp: '^Visium Human Transcriptome Probe Set v2.1.0$',
        },
        {
          key: 'Visium Mouse Transcriptome Probe Set',
          exp: '^Visium Mouse Transcriptome Probe Set$',
        },
        {
          key: 'Visium Mouse Transcriptome Probe Set v2.0',
          exp: '^Visium Mouse Transcriptome Probe Set v2.0$',
        },
        {
          key: 'Visium Mouse Transcriptome Probe Set v2.1.0',
          exp: '^Visium Mouse Transcriptome Probe Set v2.1.0$',
        },
        {
          exp: '.*',
          key: 'unknown',
        },
      ],
    },
  },

  martian_metrics: function(success_groups, failure_groups, no_perf_success_groups=[])
    local perf_groups = std.filter(function(g) !std.member(no_perf_success_groups, g), success_groups);
    [
      {
        groups: success_groups + failure_groups,
        name: 'runtime',
        special: 'runtime',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'core_hours',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'core_hours',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'in_bytes',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'in_bytes',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'out_bytes',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'out_bytes',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'in_bytes_peak',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'in_bytes_peak',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'out_bytes_peak',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'out_bytes_peak',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'total_bytes',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'total_bytes',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'output_bytes',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'output_bytes',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'vdr_bytes',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'vdr_bytes',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'forks',
            'fork_stats',
            'usertime',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'utime',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          json_path: [
            'maxbytes',
          ],
          path: [
            '_perf',
          ],
        },
        groups: perf_groups,
        name: 'disk_hwm',
        values: {
          exponential: {
            min: 1,
            multiple: 2,
          },
        },
      },
      {
        groups: perf_groups,
        name: 'perf',
        special: 'perf_stats',
      },
      {
        groups: success_groups + failure_groups,
        mro_flag: '--jobmode',
        name: 'jobmode',
        values: {
          match: [
            {
              exp: 'local',
            },
            {
              exp: 'lsf',
            },
            {
              exp: 'pbspro',
            },
            {
              exp: 'sge',
            },
            {
              exp: 'slurm',
            },
            {
              exp: 'torque',
            },
          ],
        },
      },
    ],
  named_re_matcher: function(k, insensitive=false) {
    exp: '%s%s$' % [if insensitive then '(?i)^' else '^', k],
    key: k,
  },
  prefix_re_matcher: function(k, insensitive=false) {
    exp: '%s%s' % [if insensitive then '(?i)^' else '^', k],
    key: k,
  },
}
