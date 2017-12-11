from datetime import datetime


class Logger(object):
    
    def log(self, txt): 
        '''Prints a time-stamped text string.
        
        Args:
            txt (str): text string to print.
        '''
        
        timestamped_log_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ": " + txt
        print timestamped_log_str
