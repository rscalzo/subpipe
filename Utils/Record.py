#!/usr/bin/env python

# ============================================================================
# RS 2012/02/13:  Implementation of Record class hierarchy
# ----------------------------------------------------------------------------
# This code defines a class hierarchy for constructing "record"-like objects,
# i.e., objects with attributes, from external data, including ASCII files.
# ============================================================================

import re
import os
import pyfits
import psycopg2
import psycopg2.extras


# Somewhere in here we need to define an exception for bad initialization.
# The below is an approximation but probably a crappy one.

class BadInit(Exception):
    def __init__(self, msg, **kwargs):
        self.msg = msg
        self.kwargs = kwargs
    def __str__(self):
        return self.msg


# Base class for all the APIs below
class BaseRecord(object):  pass


class KeywordRecord(BaseRecord):
    """
    The simplest API imaginable for a record-like thing.  It can be used
    to override properties read in from other APIs.
    """

    # This is a list of tuples (field, default).
    # Example: ("x", "0.0")
    # Overload as a class attribute in derived classes to implement a record.
    _keyword_fields = [ ]

    def _init_from_keywords(self, **kwargs):
        for f in self._keyword_fields:
            if kwargs.has_key(f[0]):  setattr(self,f[0],kwargs[f[0]])
            elif hasattr(self,f[0]):  pass
            else:                     setattr(self,f[0],f[1])

    def __init__(self, **kwargs):
        self._init_from_keywords(**kwargs)


class RegexRecord(BaseRecord):
    """
    Represents the results of parsing groups from a regular expression.
    """

    _regex = ""
    _regex_fields = [ ]

    def _init_from_regex(self, ascii_line):
        match = re.search(self._regex, ascii_line)
        if match == None:  return
        vals = match.groups()
        try:
            for f in self._regex_fields:  setattr(self,f[0],vals[f[1]])
        except IndexError:
            msg = "Bad index filling {0} from regular expression\n" \
                  .format(self.__class__.__name__) \
                + "field = {0}, index = {1}, regex = {2}" \
                  .format(f[0],f[1],self._regex)
            raise BadInit(msg)
        except:
            raise BadInit("Uncaught exception reading {0} from regex"
                         .format(self.__class__.__name__))

            
class AsciiRecord(BaseRecord):
    """
    Represents a line of a multi-column text file with clear delimiters.
    """

    # This is a list of tuples (field, index, informat, outformat, default)
    # governing how the record is read and written.  Basically,
    #     input:   self.field = informat(column[index]) or default
    #     output:  outstring = outformat.format(self.field)
    # Example:  [ "name", 1, (lambda x: x.strip()), '{0:s}', '' ]
    # Overload as a class attribute in derived classes to implement a record.
    # RS 2012/05/25:  Set _ascii_separator to None by default, so that items
    # in a record can be separated by arbitrary amounts of whitespace.
    _ascii_separator = None
    _ascii_fields = [ ]

    def _init_from_ascii(self, ascii_line):
        cols = ascii_line.split(self._ascii_separator)
        try:
            for f in self._ascii_fields:
                # RS 2012/06/03:  If field is empty, set it to None
                if re.search("\S+", cols[f[1]]): val = f[2](cols[f[1]])
                else:  val = None
                # FY 2014/07/09: If val is None, set field to default
                if val is None: setattr(self, f[0], f[4]) 
                else: setattr(self, f[0], val)
        except:
            raise BadInit("Badly formatted {0} ASCII table line:\n{1}"
                          .format(self.__class__.__name__, ascii_line),
                          ascii_line=ascii_line)

    def __init__ (self, ascii_line):
        self._init_from_ascii(ascii_line)

    def asline(self):
        return self._ascii_separator.join \
            ([f[3].format(getattr(self,f[0])) for f in self._ascii_fields])

    def header(self):
        return self._ascii_separator.join \
            ([re.sub(r'[a-rtz]','s',re.sub(r'\.\d+','',f[3])).format(f[0]) for f in self._ascii_fields])

    @classmethod
    def read_ascii_file(cls, fname):
        objlist = [ ]
        with open(fname) as file:
            for line in file:
                if line[0] == '#':  continue
                if len(line.strip())==0: continue
                try:
                    objlist.append(cls(line))
                except:
                    pass
        return objlist

    @staticmethod
    def write_ascii_file(fname, objlist,header=False,append=False):
        if append and os.path.exists(fname): 
            mode="a"
        else: 
            mode="w"
        with open(fname,mode) as file:
            if header and len(objlist)>0:
                file.write("{0}\n".format(objlist[0].header()))
            file.write("\n".join([obj.asline() for obj in objlist]))
            if append:
                file.write("\n")

class FITSRecord(BaseRecord):
    """
    Represents a row of a FITS binary table.
    """

    # This is a list of tuples (field, colname, type, units, format, default)
    # governing how the record is read in from the FITS file.  The fields
    # should correspond to the CFITSIO conventions for FITS binary tables.
    # Example:  [ "xsub", "X_IMAGE", "1E", "pixel", "F7.2",  0.0 ]
    # Overload as a class attribute in derived classes to implement a record.
    _fits_fields = [ ]

    def _init_from_fitsrow(self, row, colnames):
        try:
            for f in self._fits_fields:
                setattr(self,f[0],row[colnames.index(f[1])])
        except IndexError:
            msg = "Bad index filling {0} from FITS file\n" \
                  .format(self.__class__.__name__) \
                + "len(row) = {0}, ".format(len(row)) \
                + "len(colnames) = {0}, ".format(len(colnames)) \
                + "len(_fields) = {0}".format(len(self._fits_fields))
            raise BadInit(msg)
        except Exception as e:
            raise e
            raise BadInit("Uncaught exception reading {0} from FITS file"
                         .format(self.__class__.__name__))

    def __init__(self, row, colnames):
        self._init_from_fitsrow(row, colnames)

    @classmethod
    def read_fits_file(cls, fname, ext=1):
        objlist = [ ]
        with pyfits.open(fname) as ptr:
            header = ptr[0].header
            tbdata, colnames = ptr[ext].data, ptr[ext].columns.names
            if not tbdata is None:
                for r in tbdata:  objlist.append(cls(r,colnames))
        return objlist, header

    @staticmethod
    def write_fits_file(fitsname, objlist, header=None, keywords=[ ]):
        """Writes a series of FITSRecords to disk.

        RS 2012/08/22:  Now allows objlist to be an iterable of FITSRecords
        possibly of different types.  Each collection of FITSRecords of the
        same type are written to a different extension in the order in which
        they appear in objlist.
        """
        # Determine which types of objects we're dealing with.
        if len(objlist) == 0:
            return
        try:
            types = [ ]
            for s in objlist:
                if not isinstance(s, FITSRecord):  raise TypeError()
                if s.__class__ not in types:  types.append(s.__class__)
        except TypeError:
            raise BadInit("objlist must be an iterable of FITSRecords")
        # If a pyfits header object is passed in, copy it to the primary HDU.
        # Same with any one-off keywords the user would like to supply,
        # which should appear as a list of (keyword, value, comment) tuples.
        hdulist = [pyfits.PrimaryHDU(header=header)]
        for kw in keywords:  hdulist[0].header.update(*kw)
        # Now create a new table extension for each type.
        for t in types:
            col_objs = [s for s in objlist if isinstance(s, t)]
            if len(col_objs) == 0:  continue
            try:
                # This parses the data into columns for a table to be written
                # out all in one go, rather than line by line.
                cols = [pyfits.Column
                            ( name   = f[1],
                              format = f[2],
                              unit   = f[3],
                              disp   = f[4],
                              array  = [getattr(s, f[0]) for s in col_objs] )
                        for f in t._fits_fields]
            except (TypeError, IndexError):  continue
            hdulist.append(pyfits.new_table(pyfits.ColDefs(cols)))
            hdulist[-1].header.update("CLASSNAM", t.__name__,
                                       "Class name of table row objects")
        # Finally, write the primary header with extensions to disk, but only
        # if we have at least one extension so we don't end up with weird
        # corrupted files.
        if len(hdulist) > 1:
            pyfits.HDUList(hdulist).writeto(fitsname, clobber=True)
        return


class PGRecord(BaseRecord):
    """
    Represents a row of a Postgres table.  Uses psycopg.
    """

    _dbname = None
    _dbhost = None
    _dbuser = None
    _dbpw = None
    _sqlquery = ''

    # This is a list of tuples (field, dbcolname, default).
    # Example: ("gmag_err", "e_g", "0.0")
    # Basically the same as KeywordRecord, but with the possibility that
    # the field name in the database might not be what the user calls it.
    _sql_fields = [ ]

    def _init_from_sqlrow(self, **kwargs):
        for f in self._sql_fields:
            if kwargs.has_key(f[1]):  setattr(self,f[0],kwargs[f[1]])
            elif hasattr(self,f[0]):  pass
            else:                     setattr(self,f[0],f[2])

    def __init__(self, **kwargs):
        self._init_from_sqlrow(**kwargs)

    @classmethod
    def pgquery(cls, *args, **kwargs):
        conn = psycopg2.connect(host=cls._dbhost, database=cls._dbname,
                                user=cls._dbuser, password=cls._dbpw)
        cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        if 'verbose' in kwargs and kwargs['verbose']:
            print cursor.mogrify(cls._sqlquery, args)
        cursor.execute(cls._sqlquery, args)
        results = [cls(**(dict(row))) for row in cursor.fetchall()]
        # Clean up
        cursor.close()
        conn.close()
        return results
