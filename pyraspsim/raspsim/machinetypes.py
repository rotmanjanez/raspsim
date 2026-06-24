from __future__ import annotations

class iNMeta(type):
    """
    A metaclass for creating N-bit integer classes.
    """
    def __new__(mcs, name, bases, dct, n):
        dct['n'] = n
        return super().__new__(mcs, name, bases, dct)

    def __init__(cls, name, bases, dct, n):
        super().__init__(name, bases, dct)
        del n


class iN(metaclass=iNMeta, n=0):
    """
    A class representing an N-bit integer.
    """
    n: int

    def __init__(self, value: int | None = None):
        if value is None:
            value = 0

        self.value = value % (2**self.n)

    def __repr__(self):
        negrep = ""
        if self.value >> (self.n - 1):
            negrep = ", -" + str(2**self.n - self.value)

        return f"i{self.n}({self.value}, {self.value:#x}{negrep})"

    def __add__(self, other: iN | int):
        if isinstance(other, int):
            return self.__class__(self.value + other)
        return self.__class__(self.value + other.value)

    def __sub__(self, other: iN | int):
        if isinstance(other, int):
            return self.__class__(self.value - other)
        return self.__class__(self.value - other.value)

    def __mul__(self, other: iN):
        return self.__class__(self.value * other.value)

    def __truediv__(self, other: iN):
        return self.__class__(self.value // other.value)

    def __floordiv__(self, other: iN):
        return self.__class__(self.value // other.value)

    def __mod__(self, other: iN):
        return self.__class__(self.value % other.value)

    def __pow__(self, other: iN):
        return self.__class__(self.value**other.value)

    def __lshift__(self, other: iN):
        return self.__class__(self.value << other.value)

    def __rshift__(self, other: iN):
        return self.__class__(self.value >> other.value)

    def __and__(self, other: iN):
        return self.__class__(self.value & other.value)

    def __or__(self, other: iN):
        return self.__class__(self.value | other.value)

    def __xor__(self, other: iN):
        return self.__class__(self.value ^ other.value)

    def __neg__(self):
        return self.__class__(-self.value)

    def __pos__(self):
        return self.__class__(+self.value)

    def __invert__(self):
        return self.__class__(~self.value)

    def __eq__(self, other: object):
        if isinstance(other, int):
            return self.value == other
        if isinstance(other, iN):
            return self.value == other.value
        raise TypeError(f"Cannot compare i{self.n} with {type(other)}")

    def __ne__(self, other: object):
        if isinstance(other, int):
            return self.value != other
        if isinstance(other, iN):
            return self.value != other.value
        raise TypeError(f"Cannot compare i{self.n} with {type(other)}")

    def __lt__(self, other: object):
        if isinstance(other, int):
            return self.value < other
        if isinstance(other, iN):
            return self.value < other.value
        raise TypeError(f"Cannot compare i{self.n} with {type(other)}")

    def __le__(self, other: object):
        if isinstance(other, int):
            return self.value <= other
        if isinstance(other, iN):
            return self.value <= other.value
        raise TypeError(f"Cannot compare i{self.n} with {type(other)}")

    def __gt__(self, other: object):
        if isinstance(other, int):
            return self.value > other
        if isinstance(other, iN):
            return self.value > other.value
        raise TypeError(f"Cannot compare i{self.n} with {type(other)}")

    def __ge__(self, other: object):
        if isinstance(other, int):
            return self.value >= other
        if isinstance(other, iN):
            return self.value >= other.value
        raise TypeError(f"Cannot compare i{self.n} with {type(other)}")

    def __bool__(self):
        return bool(self.value)

    def __int__(self):
        return self.value


class i8(iN, metaclass=iNMeta, n=8):
    """
    A class representing an 8-bit integer.
    """


class i16(iN, metaclass=iNMeta, n=16):
    """
    A class representing a 16-bit integer.
    """


class i32(iN, metaclass=iNMeta, n=32):
    """
    A class representing a 32-bit integer.
    """


class i64(iN, metaclass=iNMeta, n=64):
    """
    A class representing a 64-bit integer.
    """
