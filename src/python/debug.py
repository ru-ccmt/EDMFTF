#!/usr/bin/env python


from types import FileType
from inspect import trace
from traceback import print_exc
from sys import stderr
import math  # log and floor


def print_stacktrace(exception=None, file=stderr, context=10, printshort=True):
    '''Prints detailed backtrace including local variables in
    each stack frame.  Should be called inside except block of
    a try-except statement.

    exception  -- exception object caught
    file       -- file object to write stacktrace
    context    -- how many lines of code context to print
    printshort -- whether or not to print short stacktrace before long one

    Example:
    try:
        do_something
    except Exception, e:
        print_stacktrace(e)
    '''
    # Handle possible errors in arguments to THIS function
    if type(file) is not FileType:
        errstr = 'WARNING: print_stacktrace called with non-file %s.  Writing to stderr.' % str(file)
        file = stderr
        print >> file, errstr
        print >> file
    elif file.closed:
        errstr = 'WARNING: print_stacktrace called with a closed file %s.  Writing to stderr.' % str(file)
        file = stderr
        print >> file, errstr
        print >> file

    if type(context) is not int:
        print >> file, 'WARNING: print_stacktrace called with context = %s which is not an integer.  Resetting to default.' % str(context)
        context = 10
    if context < 1:
        print >> file, 'WARNING: print_stacktrace called with context = %d < 1.  Resetting to 1.' % context
        context = 1

    # if desired, print short stacktrace
    if printshort:
        print >> file
        print_exc(file=file)

    # Begin processing stack
    print >> file
    print >> file, 'Detailed traceback (most recent call last):'

    framelist = trace(context=context)  # current stackframe is included
    FILL_COLUMN = 79

    try:
        for i, (frame, filename, curline, funcname, srccode, icurline) in enumerate(framelist):
            header = 'File "%s", line %d, in %s' % (filename, curline, funcname)
            print >> file, header, '-'*max((FILL_COLUMN-len(header)), 0)

            # print code block with line numbers
            caret = ['>' if i == icurline else ' ' for i in range(len(srccode))]
            linenums = [curline+(i-icurline) for i in range(len(srccode))]
            format = '%' + str(int(math.floor(math.log(max(linenums),10))+1)) + 'd '
            linestrs = [format % num for num in linenums]
            print >> file, ''.join([caret+num+line for caret,num,line in zip(caret,linestrs,srccode)])

            # print all local variables
            print >> file, 'Local variables:'
            for k,v in frame.f_locals.iteritems():
                print >> file, k, type(v), '=', v

            print >> file
    finally:
        del framelist  # avoid leakage caused by having frame references sitting around

    # if desired, print exception
    if exception is not None:
        print >> file, '%s: %s' % (exception.__class__.__name__, str(exception))

    file.flush()



if __name__ == '__main__':
    def f(a,b):
#        execfile('test.py')
        return a[0]+b

    def bad(x, y=3):
        return f(x,y)

    def main():
        try:
            x = 1
#            1/0
            print bad(3)
        except Exception, e:
            print_stacktrace(e)
    main()
