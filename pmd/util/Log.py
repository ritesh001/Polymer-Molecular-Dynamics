import logging


class CustomFormatter(logging.Formatter):

    grey = '\x1b[38;20m'
    yellow = '\x1b[33;20m'
    red = '\x1b[31;20m'
    bold_red = '\x1b[31;1m'
    reset = '\x1b[0m'
    format = '[%(levelname)s] %(message)s'

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(fmt=log_fmt,
                                      datefmt='%Y-%m-%d, %H:%M:%S %p')
        return formatter.format(record)


class Pmdlogging():
    '''Pmd's custom logging class

        Attributes:
            None
    '''

    LOGGER = logging.getLogger()
    handler = logging.StreamHandler()
    handler.setFormatter(CustomFormatter())

    LOGGER.setLevel(logging.INFO)
    LOGGER.addHandler(handler)

    @classmethod
    def info(cls, text: str):
        cls.LOGGER.info(text)

    @classmethod
    def warning(cls, text: str):
        cls.LOGGER.warning(text)
