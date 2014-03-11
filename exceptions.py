__author__ = 'kablag'


class FastRFLPError(Exception):
    def __init__(self, code):
        self.code = code

    def __str__(self):
        return repr(self.code)


class TemplateRecognitionError(FastRFLPError): pass


class GetSNPFromSequenceError(FastRFLPError): pass


class DigestError(Exception): pass


class UpdateRebaseError(FastRFLPError): pass

class SitesCollisionError(Exception): pass