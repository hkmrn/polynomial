class PolynomialError(Exception):
    """Base class for all polynomial-related exceptions"""
    pass


class PolynomialParseError(PolynomialError):
    """Raised when from_string fails to parse a valid polynomial"""
    def __init__(self, term):
        self.term = term
        super().__init__(f"Malformed polynomial term: '{term}'")


class PolynomialTypeError(PolynomialError):
    """Raised when an unsupported type is used in operations"""
    def __init__(self, op, obj):
        typename = type(obj).__name__
        super().__init__(f"Cannot {op} polynomial with object of type '{typename}'")


class PolynomialDomainError(PolynomialError):
    """Raised for invalid mathematical inputs (e.g. negative exponentiation)"""
    def __init__(self, msg="Invalid operation on polynomial"):
        super().__init__(msg)
