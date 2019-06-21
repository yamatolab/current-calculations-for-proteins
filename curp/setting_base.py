"""setting-base - Methods and class of the setting object

Classes:
SettingBase
Section
TypeBase
|+ Int
|+ Float
|+ String
|+ Bool
|+ Choice
|+ List
|+ File
"""
import os
import glob

from exception import CurpException

class SectionNotDefined(CurpException): pass
class InvalidKeyword(CurpException): pass
class InvalidType(CurpException): pass
class InvalidValue(CurpException): pass
class KeywordNotFound(CurpException): pass
class ValueCanNotConvert(CurpException): pass
class FileNotFound(CurpException): pass

################################################################################
class SettingBase:

    def __init__(self, config, check=True):

        # initialize each of sections.
        secnames = [ aname for aname in dir(self) if not aname.startswith('_') ]
        for sname in secnames:
            self[sname].set_name(sname)

        # check the given items in each section and convert to valid items.
        for secname in config.sections():
            try:
                section = getattr(self, secname)
            except AttributeError:
                raise SectionNotDefined(secname)

            pairs = config.items(secname)
            section.set_items(dict(pairs), check)


    def __getitem__(self, attrname):
        return getattr(self, attrname)

    def __str__(self):
        secnames = [ aname for aname in dir(self)
                if not aname.startswith('_') ]
        return "\n\n".join( str(self[sname]) for sname in secnames )

    def __format__(self, fmt_type="rst"):
        if fmt_type == "rst":
            secnames = [ aname for aname in dir(self)
                    if not aname.startswith('_') ]
            return "\n\n".join( format(self[sname], fmt_type)
                    for sname in secnames )
        else:
            return ""


################################################################################
class Section:

    def __init__(self, description='', **kwds):
        self.__desc = description
        self.__key_to_typeobj = kwds
        self.__name = None
        # initialize each of types
        for key, type_obj in self.__key_to_typeobj.items():
            type_obj.set_name(key)

    def __str__(self):
        seclines = "[{sname}]\n\t{desc}\n\n".format(
                sname=self.__name, desc=self.__desc)

        return seclines + "\n\n".join(str(typeobj)
                for typeobj in self.__key_to_typeobj.values())

    def __format__(self, fmt_type="rst"):
        if fmt_type == "rst":
            seclines = "{sname} section\n{secline}\n{desc}\n\n".format(
                    sname=self.__name, secline=(8+len(self.__name))*'~',
                    desc=self.__desc)

            return seclines + "\n\n".join(format(typeobj, fmt_type)
                    for typeobj in self.__key_to_typeobj.values())
        else:
            return ""

        # return '\n\n'.join( ("{sname}\n{secline}\n{desc}\n\n{kwds}"
                # .format(sname=sname, secline=len(sname)*'~',
                    # desc=kwds['description'], kwds=kwds)
                # for sname, kwds in name_to_secstrs.items() ) )


        # for key, type_obj in self.__key_to_typeobj.items():
            # lines.append( '{key}\n {value}'.format(key=key,value=type_obj) )
        # return '\n'.join(lines)

    def set_name(self, name):
        self.__name = name

    def set_items(self, kwds, check=True):
        # check whether keyword is valid or not in the config file.
        for key in kwds:
            self.check_invalid(key)

        # check whether the keyword required in definition exists or not.
        if check:
            for key in self.__key_to_typeobj:
                self.check_require(key, kwds)

        # set value
        for key, value in kwds.items():
            try:
                type_obj = self.__key_to_typeobj[key]
                type_obj.set(value)
            except FileNotFound:
                if check:
                    mes = 'key : {}, value : {}'.format(key, value)
                    raise FileNotFound(mes)
                else:
                    continue

            except InvalidValue:
                mes = 'key : {}, value : {}'.format(key, value)
                raise InvalidValue(mes)
            except ValueCanNotConvert:
                mes = 'key : {}, value : {}'.format(key, value)
                raise ValueCanNotConvert(mes)

    def check_invalid(self, key):
        """Check whether keyword is valid or not in the config file."""
        if key not in self.__key_to_typeobj:
            raise InvalidKeyword(key)

    def check_require(self, key, kwds):
        """Check whether the keyword required in definition exists or not."""
        type_obj = self.__key_to_typeobj[key]
        if type_obj.require():
            if key not in kwds:
                raise KeywordNotFound(key)

    def copy(self, **new_kwds):
        import copy
        old_kwds = copy.deepcopy(self.__key_to_typeobj)
        old_kwds.update(new_kwds)
        return Section(**old_kwds)

    def get_type_items(self):
        return self.__key_to_typeobj

    def get_type(self, attrname):
        self.check()
        return self.__key_to_typeobj[attrname]

    def __getattr__(self, attrname):
        if attrname in self.__key_to_typeobj:
            return self.__key_to_typeobj[attrname].get()
        else:
            raise AttributeError

    def __getitem__(self, attrname):
        if attrname in self.__key_to_typeobj:
            return self.__key_to_typeobj[attrname].get()
        else:
            raise IndexError


################################################################################
from abc import ABCMeta, abstractmethod, abstractproperty
class TypeBase:

    __metaclass__ = ABCMeta

    def __init__(self, default, require, desc='', parser=None):
        self.__desc = desc
        self.__require = require
        self.__parser = parser
        self.__parsed_value = None
        self.__name = None
        self.set_default(default)

    def __eq__(self, type_obj):
        if self.__class__.__name__ != type_obj.__class__.__name__:
            return False
        if self.get() == type_obj.get():
            return True
        else:
            return False

    def __str__(self):
        return "{key} = {value} # {typename}; {desc}".format(
                key=self.__name, value=self._get_value(),
                typename=self.typename, desc=self.__desc)

    def __format__(self, fmt_type="rst"):
        if fmt_type == "rst":
            if self._get_value() is "":
                val = "none"
            elif self._get_value() is True:
                val = "yes"
            elif self._get_value() is False:
                val = "no"
            else:
                val = self._get_value()

            key_val_lines = (
                    "\n**{key} = {value}** (default) : {typename}\n\n"
                    + "    {desc}").format( key=self.__name, value=val,
                            typename=self.typename, desc=self.__desc)

        else:
            key_val_lines = ""

        return key_val_lines

    def set_name(self, name):
        self.__name = name

    def set_default(self, value):
        if not self.check_default(value):
            raise InvalidValue(value)
        self.__value = self.convert_default(value)

    def set(self, value):
        if not self.check(value):
            raise ValueCanNotConvert(value)
        self.__value = self.convert(value)

    def get(self, use_parse=True):
        if use_parse:
            if self.__parser is None:
                return self.__value
            else:
                if self.__parsed_value is None:
                    self.__parsed_value = self.__parser(self.__value)
                return self.__parsed_value
        else:
            return self.__value

    def get_parser(self):
        return self.__parser

    def require(self):
        return self.__require

    def _get_value(self):
        """Get the representable value."""
        return self.get()

    @abstractproperty
    def typename(self):
        return

    @abstractmethod
    def check_default(self, value):
        return

    @abstractmethod
    def convert_default(self, value):
        return

    @classmethod
    @abstractmethod
    def check(cls, raw_value):
        """Check raw value.  This method must use only the class variables."""
        return

    @classmethod
    @abstractmethod
    def convert(cls, raw_value):
        """Return converted value.
        This method must use only the class variables.
        """
        return


################################################################################
class Int(TypeBase):

    def __init__(self, default, require=True, desc='', parser=None):
        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self): return self.__class__.__name__

    def check_default(self, value):
        return isinstance(value, int)

    def convert_default(self, value):
        return value

    @classmethod
    def check(cls, raw_value):
        try:
            int(raw_value)
            return True
        except ValueError:
            return False

    @classmethod
    def convert(cls, raw_value):
        return int(raw_value)


class Float(TypeBase):

    def __init__(self, default, require=True, desc='', parser=None):
        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self): return self.__class__.__name__

    def check_default(self, value):
        return isinstance(value, float)

    def convert_default(self, value):
        return value

    @classmethod
    def check(cls, raw_value):
        try:
            float(raw_value)
            return True
        except ValueError:
            return False

    @classmethod
    def convert(cls, raw_value):
        return float(raw_value)


class String(TypeBase):

    def __init__(self, default, require=True, desc='', parser=None):
        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self): return self.__class__.__name__

    def check_default(self, value):
        return self.check(value)

    def convert_default(self, value):
        return value.strip()

    @classmethod
    def check(cls, raw_value):
        return not raw_value.isspace()

    @classmethod
    def convert(cls, raw_value):
        return raw_value


class Bool(TypeBase):

    _trues  = ['true','True','TRUE','yes','Yes','YES','on','On','ON']
    _falses = ['false','False','FALSE','no','No','NO','off','Off','OFF']

    def __init__(self, default, require=True, desc='', parser=None):
        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self): return self.__class__.__name__

    def check_default(self, value):
        return isinstance(value, bool)

    def convert_default(self, value):
        return value

    @classmethod
    def check(cls, raw_value):
        if raw_value in cls._trues:  return True
        if raw_value in cls._falses: return True
        return False

    @classmethod
    def convert(cls, raw_value):
        if raw_value in cls._trues:  return True
        if raw_value in cls._falses: return False


class Choice(TypeBase):

    def __init__(self, default, values, value_type, require=True,
            desc='', parser=None):
        self.__Type = value_type

        # store given values
        if len(values) == 0:
            raise ValueCanNotConvert(values)
        self.__values = [ self.__Type( self.__Type.convert(raw_value) )
                          for raw_value in values ]

        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self):
        choice_str = '|'.join( [v.get() for v in self.__values] )
        return "Choice[{values}]".format(values = choice_str)

    def get(self):
        return TypeBase.get(self).get()

    def check_default(self, value):
        return self.check(value)

    def convert_default(self, value):
        return self.convert(value)

    def check(self, raw_value):
        return self.__Type( self.__Type.convert( raw_value) ) in self.__values
        # if value in [v.get() for v in self.__values]:
        #     return True
        # else:
        #     return False

    def convert(self, raw_value):
        return self.__Type( self.__Type.convert(raw_value) )


class List(TypeBase):

    def __init__(self, default, value_type, require=True, desc='', parser=None):
        self.__Type = value_type
        self.__parsed_value = None
        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self):
        return "List[{typename}]".format(typename=self.__Type.__name__)

    def get(self):
        values = [value.get() for value in TypeBase.get(self, False)]
        if self.get_parser() is None:
            return values
        else:
            if self.__parsed_value is None:
                self.__parsed_value = self.get_parser()(values)
            return self.__parsed_value

    def check_default(self, value):
        if isinstance(value, list) or isinstance(value, tuple):
            for v in value:
                try:
                    self.__Type(v)
                except InvalidValue:
                    return False
            return True
        else:
            return False

    def convert_default(self, value):
        return [ self.__Type(v) for v in value ]

    def check(self, value):
        for rv in value.split(' '):
            if not self.__Type.check(rv):
                return False
        else:
            return True

    def convert(self, values_str):
        return [ self.__Type( self.__Type.convert(raw_value) )
                 for raw_value in values_str.split(' ') if raw_value != '']

    def _get_value(self):
        return '  '.join(str(v) for v in self.get())


class File(TypeBase):

    def __init__(self, default='', require=True, desc='',
            allow_glob=True, exists=False, parser=None):
        # self._set_defaults(default)
        self.__use_glob = allow_glob
        self.__exists = exists
        TypeBase.__init__(self, default, require, desc, parser)

    @property
    def typename(self):
        return "File"

    def _set_defaults(self, default):
        if default == '':
            self.__filenames = []
        else:
            self.set(default)

    def check_default(self, value):
        if isinstance(value, str):
            return True

        elif isinstance(value, list) or isinstance(value, tuple):
            return True

        else:
            return False

    def convert_default(self, value):
        if isinstance(value, str):
            return [value]
        else:
            files = []
            for path in value:
                if self.__use_glob:
                    for file in sorted(glob.glob(path)):
                        fn = os.path.abspath(file)
                        files.append(fn)
                else:
                    files.append(path)
            return files

    def check(self, value):
        if value == '':
            return False

        for path in value.split(' '):
            if self.__use_glob:
                for fn in sorted(glob.glob(path)):
                    abs_fn = os.path.abspath(fn)
                    if self.__exists:
                        if not os.path.exists(abs_fn):
                            raise FileNotFound(fn)
            else:
                abs_fn = os.path.abspath(path)
                if self.__exists:
                    if not os.path.exists(abs_fn):
                        raise FileNotFound(path)
        return True

    def convert(self, value):
        files = []
        if isinstance(value, str):
            path = value
            if self.__use_glob:
                for fn in sorted(glob.glob(path)):
                    files.append(os.path.abspath(fn))
            else:
                files.append(os.path.abspath(path))

        else:
            for path in value.split(' '):
                if self.__use_glob:
                    for fn in sorted(glob.glob(path)):
                        files.append(os.path.abspath(fn))
                else:
                    files.append(os.path.abspath(path))
        return files

    def _get_value(self):
        return '  '.join(v for v in self.get())


################################################################################
def test_main():

    def parse(xs):
        ys = []
        for x in xs:
            ys.append(x**2)

        return ys

    generic = Section(
        int_test = Int(
            desc = '',
            parser = lambda x: x+10,
            default = 5 ),

        float_test = Float(
            desc = '',
            default = 5.5 ),

        string_test = String(
            desc = '',
            default = 'hoge' ),

        format = Choice(
            desc = '',
            values = ['presto', 'amber'],
            value_type = String,
            default = 'presto' ),

        crd = List(
            desc = '',
            value_type = Int,
            default = [1,2,3] ),

        crd_power = List(
            desc = '',
            value_type = Int,
            parser = parse,
            default = [1,2,3] ),

        target = Choice(
            desc = '',
            values = ['trajectory', 'restart'],
            value_type = String,
            default = 'restart' ),

        bool_test = Bool(
            desc = '',
            default = True),

        file_test = File(exists=False,
            desc = '',
            default = 'flux.dat')

        )

    # print(str(generic))
    # print(generic['format'])

    generic.set_items(dict(float_test=0, file_test='run.ini *g.py'),check=False)
    print(str(generic))


    # print(generic.format)
    # print(generic.int_test)
    # print(generic.float_test)
    # print(generic.bool_test)
    # print(generic.string_test)
    # print(generic.file_test)

if __name__ == '__main__':
    test_main()
