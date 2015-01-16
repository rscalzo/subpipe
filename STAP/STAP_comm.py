import threading
import subprocess


def print_cmd_line(name, *args, **kwargs):
    print "To reproduce this step, run:"
    print "    $SUBPIPEHOME/STAP/{0}".format(name),
    for arg in args:
        if type(arg) is list:
            print " ".join(arg),
        else:
            print arg,
    # We use some parsing conventions below for keyword arguments.
    # If the value is None, we assume it has been omitted and don't print it.
    # If the value is of Boolean type, we assume it's action="store_true".
    for (kw,kwval) in kwargs.items():
        if kwval is not None:
            if type(kwval) is bool:
                if kwval:  print "--{0}".format(kw),
            elif type(kwval) is list:
                print "--{0} {1}".format(kw, " ".join(kwval)),
            elif type(kwval) is str:
                print "--{0} \'{1}\'".format(kw, kwval),
            else:
                print "--{0}={1}".format(kw, kwval),
    print; print


def psmemchk(cmdre, verbose=False):
    # RS 2012/07/04:  Use ps to check which processes are running and
    # how much memory they're taking up.  It isn't the same as profiling
    # but it'll answer questions about whether problems accumulate when
    # processes die or whether they're in there from the start.
    memtot = 0
    if verbose:  print "Subprocesses in play:"
    proc = subprocess.Popen("ps -u rscalzo -o pid,ppid,pmem,rss,cmd "
                            "| grep '{0}'".format(cmdre), shell=True,
                            stdout=subprocess.PIPE)
    stdoutstr = proc.communicate()[0]
    for line in stdoutstr.split("\n"):
        if verbose:  print "   ", line
        try:
            cols = line.split()
            memtot += float(cols[3])/1024.0
        except:
            continue
    # NB:  "pmem" below is hard-wired for maipenrai (112 GB RAM)
    pmemtot = memtot*100/114688.0
    if verbose:
        print "Using {0:.1f} MB ({1:.1f}%) of available memory".format(
               memtot, pmemtot)
    return memtot, pmemtot



def STAP_callexternal(cmd, timeout=None, combinestderr=False,
                      getstdout=False, getstderr=False, input=None,
                      stdout=None, stderr=None, verbose=True):
    """
    Call an external command with a given timeout.  Blatantly rips off some
    code in class Command (see below) from: https://gist.github.com/1306188
    Might want to consider making that code raise TimeOut if it times out.
    """

    myin, myout, myerr = None, subprocess.PIPE, subprocess.PIPE
    if stdout:  myout = stdout
    if stderr:  myerr = stderr
    if combinestderr:  myerr = subprocess.STDOUT
    if input:  myin = subprocess.PIPE

    if verbose:  print "STAP_callexternal running:  {0}".format(cmd)
    cmd = Command(cmd)
    result = cmd.run(timeout=timeout, input=input,
                     stdout=myout, stderr=myerr, stdin=myin)
    if not getstdout and not getstderr:
        return result
    elif getstdout and not getstderr:
        return result, cmd.stdoutstr
    elif getstderr and not getstdout:
        return result, cmd.stderrstr
    else:
        return result, cmd.stdoutstr, cmd.stderrstr


class Command(object):
    '''
    Enables to run subprocess commands in a different thread
    with TIMEOUT option!

    Based on jcollado's solution:
    http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
    '''

    def __init__(self, cmd):
        self.cmd = cmd.split(" ")
        self.process = None

    def run(self, timeout=None, input=None, **kwargs):
        def target(**kwargs):
            self.process = subprocess.Popen(self.cmd, **kwargs)
            self.stdoutstr, self.stderrstr = self.process.communicate(input)

        thread = threading.Thread(target=target, kwargs=kwargs)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()

        # RS 2012/06/19:  Status codes should be *signed* integers.
        status = self.process.returncode
        if status >= 128: status -= 256
        return status
