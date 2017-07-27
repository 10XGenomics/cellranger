var LZString = require('./external/lz-string-1.4.4.min.js');

var stdin = process.stdin;
var stdout = process.stdout;
var inputChunks = [];

stdin.setEncoding('utf8');

stdin.on('data', function (chunk) {
  inputChunks.push(chunk);
});

stdin.on('end', function () {
  var inputString = inputChunks.join("");
  var uncompressed = LZString.decompressFromEncodedURIComponent(inputString);
  stdout.write(uncompressed);
});
