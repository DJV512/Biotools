from functools import wraps
from inspect import signature

def validate_dna(allow_n=False):
    '''
    Decorator to inspect DNA sequences and insure validity before operating on them.
    Works on functions that accept one or multiple dna sequences, as long as their parameter names start with "dna."
    Accepts its own parameter, "allow_n" which defaults to False, but can be set to True if a function is allowed to work with DNA that contains ambiguous "N" bases.
    '''
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            sig = signature(func)
            bound = sig.bind(*args, **kwargs)

            if allow_n:
                allowed_bases = {"A", "C", "G", "T", "N"}
            else:
                allowed_bases = {"A", "C", "G", "T"}

            for name, value in bound.arguments.items():
                if name.startswith("dna"):
                    bases = set(value.upper())
                    if not bases.issubset(allowed_bases):
                        raise ValueError(f"{value} is not valid a DNA sequence.")
            
            return func(*args, **kwargs)
        return wrapper
    return decorator

