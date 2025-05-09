from exceptions import PolynomialParseError, PolynomialTypeError, PolynomialDomainError
from typing import Union

Scalar = Union[int, float]

class Polynomial:
    """
    Represents a univariate polynomial with real coefficients.
    """

    def __init__(self, *args: Scalar, var: str = "x") -> None:
        """
        Initialises a polynomial from its coefficients.

        Example: Polynomial(1, 0, -3) -> 1 - 3x^2
        Terms are in ascending order of powers.
        """
        self._coefficients = self._clean(self._strip_leading_zeros(args))
        self._var = var
        if not self._is_valid_coefficients(self.coefficients):
            raise ValueError("Polynomial coefficients must be real numbers")

    @property
    def coefficients(self) -> tuple[Scalar, ...]:
        """
        Returns the tuple of coefficients.
        """
        return self._coefficients

    @property
    def var(self) -> str:
        """
        Returns the variable name used in the string representation.
        """
        return self._var

    @var.setter
    def var(self, new: str) -> None:
        """
        Sets the variable name. Must be a string.
        """
        if not isinstance(new, str):
            raise TypeError("Polynomial variable must be a string")
        self._var = new

    @property
    def degree(self) -> int:
        """
        Returns the degree of the polynomial.
        """
        return len(self) - 1

    def __str__(self) -> str:
        """
        Returns a human-readable algebraic representation.
        """
        result = ""
        nonzero_found = False
        for i in range(len(self.coefficients)):
            if self.coefficients[i] == 0:
                continue
            if not nonzero_found:
                plus = "" if self.coefficients[i] >= 0 else "-"
            elif self.coefficients[i] >= 0:
                plus = " + "
            else:
                plus = " - "
            coefficient = abs(self.coefficients[i]) if abs(self.coefficients[i]) != 1 or i == 0 else ""
            nonzero_found = True
            if i == 0:
                result += f"{plus}{coefficient}"
            elif i == 1:
                result += f"{plus}{coefficient}{self.var}"
            else:
                result += f"{plus}{coefficient}{self.var}^{i}"
        return result or "0"

    def __repr__(self) -> str:
        """
        Returns the code-style representation of the polynomial.
        """
        return f"Polynomial({', '.join([str(i) for i in self.coefficients])}, var=\"{self.var}\")"

    def __len__(self) -> int:
        """
        Returns the number of terms (degree + 1).
        """
        return len(self.coefficients)

    def __eq__(self, other: object) -> bool:
        """
        Checks coefficient-wise equality with scalars or other polynomials.
        """
        if Polynomial._is_numeric(other):
            return self.coefficients[0] == other
        if isinstance(other, Polynomial):
            return self.coefficients == other.coefficients
        return False

    def __hash__(self) -> int:
        """
        Allows use of Polynomial as dict keys or set elements.
        """
        return hash(tuple(self.coefficients))

    def __add__(self, other: Union[Scalar, "Polynomial"]) -> "Polynomial":
        """
        Adds a scalar or another Polynomial.
        """
        if Polynomial._is_numeric(other):
            return Polynomial(self.coefficients[0] + other, *self.coefficients[1:])
        if isinstance(other, Polynomial):
            shorter, longer = sorted([self, other], key=len)
            diff = abs(len(self) - len(other))
            shorter_coefficients = shorter.coefficients + [0] * diff
            return Polynomial(*[co1 + co2 for co1, co2 in zip(shorter_coefficients, longer.coefficients)], var=self.var)
        raise PolynomialTypeError("add", other)

    def __radd__(self, other: Scalar) -> "Polynomial":
        """
        Enables scalar + Polynomial (reversed addition).
        """
        return self + other

    def __sub__(self, other: Union[Scalar, "Polynomial"]) -> "Polynomial":
        """
        Subtracts a scalar or another Polynomial.
        """
        return self + -other

    def __neg__(self) -> "Polynomial":
        """
        Returns the negated polynomial (-self).
        """
        return Polynomial(*[-x for x in self.coefficients], var=self.var)

    def __mul__(self, other: Union[Scalar, "Polynomial"]) -> "Polynomial":
        """
        Multiplies by a scalar or another polynomial.
        """
        if Polynomial._is_numeric(other):
            return Polynomial(*[co * other for co in self.coefficients], var=self.var)
        if isinstance(other, Polynomial):
            new_coefficients = [0] * (len(self) + len(other) - 1)
            for i in range(len(self.coefficients)):
                for j in range(len(other.coefficients)):
                    new_coefficients[i + j] += self.coefficients[i] * other.coefficients[j]
            return Polynomial(*new_coefficients, var=self.var)
        raise PolynomialTypeError("multiply", other)

    def __rmul__(self, other: Scalar) -> "Polynomial":
        """
        Enables scalar * Polynomial (reversed multiplication).
        """
        return self * other

    def __truediv__(self, other: Scalar) -> "Polynomial":
        """
        Divides by a scalar.
        """
        if Polynomial._is_numeric(other):
            return Polynomial(*[co / other for co in self.coefficients], var=self.var)
        raise PolynomialTypeError("divide", other)

    def __pow__(self, n: int) -> "Polynomial":
        """
        Raises the polynomial to a non-negative integer power.
        """
        if not isinstance(n, int) or n < 0:
            raise PolynomialDomainError("Cannot raise polynomial to negative power")
        result = Polynomial(1, var=self.var)
        for _ in range(n):
            result *= self
        return result

    def evaluate(self, n: Scalar) -> Scalar:
        """
        Evaluates the polynomial at a number n.
        """
        return sum(self.coefficients[i] * n ** i for i in range(len(self)))

    def copy(self) -> "Polynomial":
        """
        Returns a deep copy of the polynomial.
        """
        return Polynomial(*self.coefficients, var=self.var)

    def derivative(self) -> "Polynomial":
        """
        Returns the first derivative of the polynomial.
        """
        if len(self.coefficients) == 1:
            return Polynomial(0)
        differentiated = [i * co for i, co in enumerate(self.coefficients)][1:]
        return Polynomial(*differentiated, var=self.var)

    def integral(self, constant: Scalar = 0) -> "Polynomial":
        """
        Returns the indefinite integral with an optional constant.
        """
        integrated = [constant] + [self.coefficients[i] / (i + 1) for i in range(len(self.coefficients))]
        return Polynomial(*integrated, var=self.var)

    def compose(self, other: "Polynomial") -> "Polynomial":
        """
        Returns the composition P(Q(x)) where self is P and other is Q.
        """
        if not isinstance(other, Polynomial):
            raise PolynomialTypeError("compose", other)
        result = Polynomial(0, var=self.var)
        for i in range(len(self.coefficients)):
            result += self.coefficients[i] * other**i
        return result

    @staticmethod
    def from_roots(*args: Scalar, var: str = "x") -> "Polynomial":
        """
        Constructs a polynomial from its roots.
        """
        result = Polynomial(1, var=var)
        for root in args:
            result *= Polynomial(-root, 1, var=var)
        return result

    @staticmethod
    def from_string(string: str, var: str = "x") -> "Polynomial":
        """
        Parses a string like "x^2 - 3x + 2" into a Polynomial object.
        """
        if string.startswith("-"):
            string = "0" + string
        string = string.replace(" ", "").replace("+-", "_").replace("-", "+-").replace("_", "+-")
        terms = string.split("+")
        coefficients = {}
        for term in terms:
            if Polynomial._is_valid_term(term, var):
                coefficient, exponent = Polynomial._is_valid_term(term, var)
                if coefficients.get(exponent):
                    coefficients[exponent] += coefficient
                else:
                    coefficients[exponent] = coefficient
            else:
                raise PolynomialParseError(term)
        degree = max(coefficients)
        final_coefficients = [0] * (degree + 1)
        for i in range(degree + 1):
            if coefficients.get(i):
                final_coefficients[i] = coefficients[i]
        return Polynomial(*final_coefficients, var=var)

    @staticmethod
    def _is_valid_term(term: str, var: str) -> Union[tuple[Scalar, int], bool]:
        """
        Parses a term like '3x^2' into (coefficient, exponent), or returns False if invalid.
        """
        try:
            single = float(term)
            return single, 0
        except:
            pass
        check = term.split(f"{var}^")
        if check[0] == "":
            check[0] = 1
        elif check[0] == "-":
            check[0] = -1
        if len(check) != 2:
            if "^" not in term:
                try:
                    sign = 1
                    if term[0] == "-":
                        sign = -1
                        term = term[1:]
                        if len(term) == 1:
                            return -1, 1
                    coefficient = float(term[:-1])
                    return sign * coefficient, 1
                except:
                    pass
            return False
        try:
            coefficient = float(check[0])
            exponent = int(check[1])
        except:
            return False
        return coefficient, exponent
