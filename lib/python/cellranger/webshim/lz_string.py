# Based on the LZ-String javascript found here:
#   http://pieroxy.net/blog/pages/lz-string/index.html
#   version 1.4.4
import os
import random
import subprocess
import unittest
import tenkit.log_subprocess as tk_subproc

keyStrUriSafe = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-$"

def compressToEncodedURIComponent(input):
    if input is None:
        return ""

    return _compress(input, 6, lambda a: keyStrUriSafe[a])

def _compress(uncompressed, bits_per_char, get_char_from_int):
    if uncompressed is None:
        return ""

    context_dictionary = {}
    context_dictionaryToCreate= {}
    context_wc = ""
    context_w = ""
    context_enlargeIn = 2  # Compensate for the first entry which should not count
    context_dictSize = 3
    context_numBits = 2
    context_data = []
    context_data_val = 0
    context_data_position = 0

    def get_char_iter(s, chunk_size=1024):
        if hasattr(s, '__len__'):
            for c in s:
                yield c
        elif hasattr(s, 'read'):
            while True:
                chunk = s.read(chunk_size)
                if len(chunk) == 0:
                    break
                for c in chunk:
                    yield c
        else:
            raise Exception("Don't know how to compress object of type %s" % type(s))

    for context_c in get_char_iter(uncompressed):
        if context_c not in context_dictionary:
            context_dictionary[context_c] = context_dictSize
            context_dictSize += 1
            context_dictionaryToCreate[context_c] = True

        context_wc = context_w + context_c;
        if context_wc in context_dictionary:
            context_w = context_wc
        else:
            if context_w in context_dictionaryToCreate:
                if ord(context_w[0]) < 256:
                    for i in xrange(context_numBits):
                        context_data_val <<= 1

                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0
                            context_data.append(get_char_from_int(context_data_val))
                            context_data_val = 0
                        else:
                            context_data_position += 1

                    value = ord(context_w[0])
                    for i in xrange(8):
                        context_data_val = (context_data_val << 1) | (value & 1);
                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0
                            context_data.append(get_char_from_int(context_data_val))
                            context_data_val = 0
                        else:
                            context_data_position += 1

                        value >>= 1
                else:
                    value = 1;
                    for i in xrange(context_numBits):
                        context_data_val = (context_data_val << 1) | value;
                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0;
                            context_data.append(get_char_from_int(context_data_val))
                            context_data_val = 0;
                        else:
                            context_data_position += 1;

                        value = 0

                    value = ord(context_w[0])
                    for i in xrange(16):
                        context_data_val = (context_data_val << 1) | (value & 1)
                        if context_data_position == bits_per_char - 1:
                            context_data_position = 0;
                            context_data.append(get_char_from_int(context_data_val))
                            context_data_val = 0;
                        else:
                            context_data_position += 1

                        value >>= 1

                context_enlargeIn -= 1
                if context_enlargeIn == 0:
                    context_enlargeIn = 2**context_numBits
                    context_numBits += 1
                context_dictionaryToCreate.pop(context_w, None)
            else:
                value = context_dictionary[context_w]

                for i in xrange(context_numBits):
                    context_data_val = (context_data_val << 1) | (value & 1)
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0;
                        context_data.append(get_char_from_int(context_data_val))
                        context_data_val = 0;
                    else:
                        context_data_position += 1

                    value >>= 1

            context_enlargeIn -= 1
            if context_enlargeIn == 0:
                context_enlargeIn = 2**context_numBits
                context_numBits += 1

            # Add wc to the dictionary.
            context_dictionary[context_wc] = context_dictSize
            context_dictSize += 1
            context_w = context_c

    # Output the code for w.
    if context_w != '':
        if context_w in context_dictionaryToCreate:
            if ord(context_w[0]) < 256:
                for i in xrange(context_numBits):
                    context_data_val <<= 1
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0;
                        context_data.append(get_char_from_int(context_data_val))
                        context_data_val = 0;
                    else:
                        context_data_position += 1

                value = ord(context_w[0])
                for i in xrange(8):
                    context_data_val = (context_data_val << 1) | (value & 1);
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0
                        context_data.append(get_char_from_int(context_data_val))
                        context_data_val = 0
                    else:
                        context_data_position += 1

                    value >>= 1
            else:
                value = 1;
                for i in xrange(context_numBits):
                    context_data_val = (context_data_val << 1) | value;
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0;
                        context_data.append(get_char_from_int(context_data_val))
                        context_data_val = 0;
                    else:
                        context_data_position += 1;

                    value = 0

                value = ord(context_w[0])
                for i in xrange(16):
                    context_data_val = (context_data_val << 1) | (value & 1)
                    if context_data_position == bits_per_char - 1:
                        context_data_position = 0;
                        context_data.append(get_char_from_int(context_data_val))
                        context_data_val = 0;
                    else:
                        context_data_position += 1

                    value >>= 1

            context_enlargeIn -= 1
            if context_enlargeIn == 0:
                context_enlargeIn = 2**context_numBits
                context_numBits += 1
            context_dictionaryToCreate.pop(context_w, None)
        else:
            value = context_dictionary[context_w]

            for i in xrange(context_numBits):
                context_data_val = (context_data_val << 1) | (value & 1)
                if context_data_position == bits_per_char - 1:
                    context_data_position = 0;
                    context_data.append(get_char_from_int(context_data_val))
                    context_data_val = 0;
                else:
                    context_data_position += 1

                value >>= 1

        context_enlargeIn -= 1
        if context_enlargeIn == 0:
            context_enlargeIn = 2**context_numBits
            context_numBits += 1

    # Mark the end of the stream
    value = 2
    for i in xrange(context_numBits):
        context_data_val = (context_data_val << 1) | (value & 1);
        if context_data_position == bits_per_char - 1:
            context_data_position = 0
            context_data.append(get_char_from_int(context_data_val))
            context_data_val = 0
        else:
            context_data_position += 1

        value >>= 1

    # Flush the last char
    while True:
        context_data_val <<= 1;
        if context_data_position == bits_per_char - 1:
            context_data.append(get_char_from_int(context_data_val))
            break
        else:
            context_data_position += 1;

    return ''.join(context_data)

class LZStringUnitTest(unittest.TestCase):
    def test_hello_world(self):
        compressed = compressToEncodedURIComponent('Hello world!')
        self.assertFalse(compressed == 'Hello world!')
        decompressed = self.decompress(compressed)
        self.assertTrue(decompressed == 'Hello world!')

    def test_repeated_string(self):
        test_string = 'aaaaabaaaaacaaaaadaaaaaeaaaaa'
        compressed = compressToEncodedURIComponent(test_string)
        self.assertFalse(compressed == test_string)
        self.assertTrue(len(compressed) < len(test_string))
        decompressed = self.decompress(compressed)
        self.assertTrue(decompressed == test_string)

    def test_long_string(self):
        test_string = '';
        for i in xrange(1000):
            test_string += str(random.random()) + " "

        compressed = compressToEncodedURIComponent(test_string);
        self.assertFalse(compressed == test_string)
        self.assertTrue("=" not in compressed)
        self.assertTrue("/" not in compressed)
        self.assertTrue(len(compressed) < len(test_string))
        decompressed = self.decompress(compressed)
        self.assertTrue(decompressed == test_string)

    def test_paragraph(self):
        decompressed = 'During tattooing, ink is injected into the skin, initiating an immune response, and cells called "macrophages" move into the area and "eat up" the ink. The macrophages carry some of the ink to the body\'s lymph nodes, but some that are filled with ink stay put, embedded in the skin. That\'s what makes the tattoo visible under the skin. Dalhousie Uiversity\'s Alec Falkenham is developing a topical cream that works by targeting the macrophages that have remained at the site of the tattoo. New macrophages move in to consume the previously pigment-filled macrophages and then migrate to the lymph nodes, eventually taking all the dye with them. "When comparing it to laser-based tattoo removal, in which you see the burns, the scarring, the blisters, in this case, we\'ve designed a drug that doesn\'t really have much off-target effect," he said. "We\'re not targeting any of the normal skin cells, so you won\'t see a lot of inflammation. In fact, based on the process that we\'re actually using, we don\'t think there will be any inflammation at all and it would actually be anti-inflammatory.'
        compressed = "CIVwTglgdg5gBAFwIYIQezdGAaO0DWeAznlAFYCmAxghQCanqIAWFcR 0u0ECEKWOEih4AtqJBQ2YCkQAOaKEQq5hDKhQA2mklSTb6cAESikVMGjnMkMWUbii0ANzbQmCVkJlIhUBkYoUOBA5ew9XKHwAOjgAFU9Tc0trW10kMDAAT3Y0UTY0ADMWCMJ3TwAjNDpMgHISTUzRKzgoKtlccpAEHLyWIPS2AogDBgB3XmZSQiJkbLku3ApRcvo6Q2hi9k4oGPiUOrhR627TfFlN5FQMOCcIIghyzTZJNbBNjmgY4H1mNBB7tgAVQgLjA9wQtRIAEEnlQ4AAxfRnKDWUTEOBrFyaSyCHzoOQQPSaODmQJojxBUZoMD4EjlbLIMC2PiwTaJCxWGznCndawuOAyUzQQxBcLsXj5Ipiy7oNAxAByFFGDjMHJS50c-I2TCoiiIIF6YrkMlufyIDTgBJgeSgCAAtEMRiqkpzUr4GOERKIIDAwCg2GU2A0mpNWmsiIsXLaQPoLchtvBY5tqmxxh5iqIYkYAOqsES6prpQS8RBoOCaJDKMB28qVwwy66C5z6bgiI6EyaZP7sCgBirgJS4MVEPQZLBDiqaO60MGtlh3El13CjCg1fnhn1SBg OhgEDwHkYtCyKA1brebTZPlsCRUSaFAp2xnMuAUAoFagIbD2TxEJAQOgs2zVcZBaNBumfCgWUTKBskKTZWjAUxiQ fMtB0XAiDLLsQEORQzx7NgfGxbp4OgAoK3EARFBiABJEQCjML84FrZQGEUTZjTQDQiBIQ8VxqUCmJjS9gnuWBlzYOh8Ig5gCGKUDxm0FiiNg0gKKQKi A4-plLUPBuipEBNG3GgRItFZfD4O1yMo0x0CyKIgAAA$$"
        decomp2 = self.decompress(compressed)
        self.assertTrue(decompressed == decomp2)

    def decompress(self, compressed):
        cwd = os.path.dirname(os.path.realpath(__file__))
        p = tk_subproc.Popen(['node', 'decompress.js'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=cwd)
        stdout, _ = p.communicate(input=compressed)
        self.assertTrue(p.returncode == 0)
        return stdout

# NOTE: nodejs is required for unit tests
if __name__ == '__main__':
    unittest.main()
