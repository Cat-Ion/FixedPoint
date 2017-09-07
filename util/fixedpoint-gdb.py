import re
class MultiwordIntegerPrettyPrinter:
    def __init__(self, val):
        self.val = val
        self.n = val.type.template_argument(0)
        self.s = val['s'][0].type.sizeof

    def to_string(self):
        return str(self.to_number())

    def to_number(self):
        ss = 8*self.s
        v = [ int(self.val['s'][i]) for i in range(self.n) ]
        if v[-1] >= 2**(ss-1):
            v[-1] -= 2**ss
        return sum([ v[i] * 2**(i*ss) for i in range(self.n) ])

    def children(self):
        return [ ('numWords', self.n), ('storageSize', self.s), ('s', self.val['s']) ]

class FixedPointPrettyPrinter:
    def __init__(self, val):
        self.val = val
    
    def to_string(self):
        return str(self.to_number())

    def to_number(self):
        v = MultiwordIntegerPrettyPrinter(self.val['v']).to_number()
        e = self.val.type.template_argument(1)
        return v * 2.**(-e)

    def children(self):
        return [('integerWidth', self.val.type.template_argument(0)),
                ('fractionalWidth', self.val.type.template_argument(1)),
                ('v', self.val['v']) ]

def lookup_function(val):
   lookup_tag = val.type.strip_typedefs().tag
   if lookup_tag == None:
       return None
   if re.compile("^MultiwordInteger<.*>$").match(lookup_tag):
        return MultiwordIntegerPrettyPrinter(val)
   if re.compile("^FixedPoint<.*>$").match(lookup_tag):
        return FixedPointPrettyPrinter(val)
   return None

def register_printers(objfile):
   objfile.pretty_printers.append(lookup_function)

register_printers(gdb.current_objfile())
