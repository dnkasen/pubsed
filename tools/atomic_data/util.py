import sys
import io
import numpy as np

# a few utility functions
def peek_line(f):
    pos = f.tell()
    line = f.readline()
    f.seek(pos)
    return line


def skip_until_found(search_string, f):
    line = peek_line(f)
    while True:
        if line.find(search_string) > 0:
            break
        line = f.readline()


def to_array(*args):
    return [np.array(arg) for arg in args]


# this is for python 2 / 3 compatibility
if sys.version_info[0] == 3:
    def iteritems(d, **kw):
        return iter(d.items(**kw))
else:
    def iteritems(d, **kw):
        return iter(d.iteritems(**kw))


# this is to deal with the fact that the text files in the
# cmfgen database have different encodings
encodings = ['utf-8', 'iso-8859-1']

def open_text_file(fn, mode):
    for e in encodings:
        try:
            f = io.open(fn, mode, encoding=e)
            f.readline()
            f.seek(0)
        except UnicodeDecodeError:
            pass
        else:
            return f
