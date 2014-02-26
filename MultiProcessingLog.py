"""Python 2 logger for multiprocessing taken from
http://stackoverflow.com/questions/641420/how-should-i-log-while-using-multiprocessing-in-python
"""

import threading, logging, sys, traceback

class QueueHandler(logging.Handler):
    def __init__(self, queue):
        logging.Handler.__init__(self)
        self.queue = queue

    def send(self, s):
        self.queue.put_nowait(s)

    def _format_record(self, record):
        # ensure that exc_info and args
        # have been stringified.  Removes any chance of
        # unpickleable things inside and possibly reduces
        # message size sent over the pipe
        if record.args:
            record.msg = record.msg % record.args
            record.args = None
        if record.exc_info:
            dummy = self.format(record)
            record.exc_info = None

        return record

    def emit(self, record):
        try:
            s = self._format_record(record)
            self.send(s)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

    def close(self):
        logging.Handler.close(self)

class QueueListener(threading.Thread):
    def __init__(self, queue):
        super(QueueListener, self).__init__()
        self.queue = queue
        self.daemon = True
        self.start()

    def run(self):
        while True:
            try:
                record = self.queue.get()
                logger = logging.getLogger(record.name)
                logger.callHandlers(record)
            except (KeyboardInterrupt, SystemExit):
                raise
            except EOFError:
                break
            except:
                traceback.print_exc(file=sys.stderr)
