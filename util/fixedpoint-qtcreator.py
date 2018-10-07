from dumper import *

def qdump__FixedPoint(d, value):
    d.putNumChild(3)
    raw = [ value["v"]["s"][i].integer() for i in range( value["v"]["numWords"].integer() ) ]
    ss = value["v"]["storageSize"].integer()
    exp = [raw[i] * 2**(i * ss) for i in range(len(raw)) ]
    if raw[-1] >= 2**(ss-1):
        exp += [ -2**(ss * len(raw)) ]
    d.putValue(sum(exp) * 2**-value["fractionalWidth"].integer())
    if d.isExpanded():
        with Children(d):
            d.putSubItem("fractionalWidth", value["fractionalWidth"])
            d.putSubItem("integerWidth", value["integerWidth"])
            d.putSubItem("v", value["v"])

def qdump__MultiwordInteger(d, value):
    d.putNumChild(3)
    raw = [ value["s"][i].integer() for i in range( value["numWords"].integer() ) ]
    exp = [ raw[i] * 2**(i * value["storageSize"].integer()) for i in range(len(raw)) ]
    d.putValue(sum(exp))
    if d.isExpanded():
        with Children(d):
            d.putSubItem("numWords", value["numWords"])
            d.putSubItem("storageSize", value["storageSize"])
            d.putSubItem("s", value["s"])

