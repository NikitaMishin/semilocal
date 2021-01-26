class Logger:
    def __init__(self, log_file):
        self.file = open(log_file, 'wt', 1)

    def write(self, msg):
        # print(msg)
        self.file.write(msg + '\n')

    def close(self):
        self.file.close()
