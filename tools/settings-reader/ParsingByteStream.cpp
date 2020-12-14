//C-  -*- C++ -*-
//C- -------------------------------------------------------------------
//C- DjVuLibre-3.5
//C- Copyright (c) 2002  Leon Bottou and Yann Le Cun.
//C- Copyright (c) 2001  AT&T
//C-
//C- This software is subject to, and may be distributed under, the
//C- GNU General Public License, either Version 2 of the license,
//C- or (at your option) any later version. The license should have
//C- accompanied the software or you may obtain a copy of the license
//C- from the Free Software Foundation at http://www.fsf.org .
//C-
//C- This program is distributed in the hope that it will be useful,
//C- but WITHOUT ANY WARRANTY; without even the implied warranty of
//C- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//C- GNU General Public License for more details.
//C- 
//C- DjVuLibre-3.5 is derived from the DjVu(r) Reference Library from
//C- Lizardtech Software.  Lizardtech Software has authorized us to
//C- replace the original DjVu(r) Reference Library notice by the following
//C- text (see doc/lizard2002.djvu and doc/lizardtech2007.djvu):
//C-
//C-  ------------------------------------------------------------------
//C- | DjVu (r) Reference Library (v. 3.5)
//C- | Copyright (c) 1999-2001 LizardTech, Inc. All Rights Reserved.
//C- | The DjVu Reference Library is protected by U.S. Pat. No.
//C- | 6,058,214 and patents pending.
//C- |
//C- | This software is subject to, and may be distributed under, the
//C- | GNU General Public License, either Version 2 of the license,
//C- | or (at your option) any later version. The license should have
//C- | accompanied the software or you may obtain a copy of the license
//C- | from the Free Software Foundation at http://www.fsf.org .
//C- |
//C- | The computer code originally released by LizardTech under this
//C- | license and unmodified by other parties is deemed "the LIZARDTECH
//C- | ORIGINAL CODE."  Subject to any third party intellectual property
//C- | claims, LizardTech grants recipient a worldwide, royalty-free, 
//C- | non-exclusive license to make, use, sell, or otherwise dispose of 
//C- | the LIZARDTECH ORIGINAL CODE or of programs derived from the 
//C- | LIZARDTECH ORIGINAL CODE in compliance with the terms of the GNU 
//C- | General Public License.   This grant only confers the right to 
//C- | infringe patent claims underlying the LIZARDTECH ORIGINAL CODE to 
//C- | the extent such infringement is reasonably necessary to enable 
//C- | recipient to make, have made, practice, sell, or otherwise dispose 
//C- | of the LIZARDTECH ORIGINAL CODE (or portions thereof) and not to 
//C- | any greater extent that may be necessary to utilize further 
//C- | modifications or combinations.
//C- |
//C- | The LIZARDTECH ORIGINAL CODE is provided "AS IS" WITHOUT WARRANTY
//C- | OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
//C- | TO ANY WARRANTY OF NON-INFRINGEMENT, OR ANY IMPLIED WARRANTY OF
//C- | MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
//C- +------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#if NEED_GNUG_PRAGMAS
# pragma implementation
#endif

#include "ParsingBytestream.h"

//#include "IW44Image.h"
#include "GOS.h"
#include "GString.h"
//#include "DjVuDocEditor.h"
//#include "DjVuDumpHelper.h"
//#include "DjVuMessageLite.h"
//#include "BSByteStream.h"
//#include "DjVuText.h"
//#include "DjVuAnno.h"
//#include "DjVuInfo.h"
//#include "IFFByteStream.h"
//#include "DataPool.h"
//#include "DjVuPort.h"
//#include "DjVuFile.h"
//#include "DjVmNav.h"
//#include "common.h"

//static bool modified = false;
//static bool verbose = false;
//static bool save = false;
//static bool nosave = false;
//static bool utf8 = false;

static unsigned char utf8bom[] = { 0xef, 0xbb, 0xbf };

//static GUTF8String
//ToNative(GUTF8String s)
//{
//  if (utf8)
//    return s;
//  // fake the damn GUTF8/GNative type check
//  GNativeString n = s;
//  return GUTF8String((const char*)n);
//}



// --------------------------------------------------
// PARSING BYTESTREAM
// --------------------------------------------------

ParsingByteStream::ParsingByteStream(const GP<ByteStream> &xgbs)
  : gbs(xgbs),bs(*gbs), bufpos(1), bufend(1), goteof(false)
{ 
}

int 
ParsingByteStream::eof() // aka. feof
{
  if (bufpos < bufend) 
    return false;
  if (goteof)
    return true;
  bufend = bufpos = 1;
  while (bs.read(buffer+bufend,1) && ++bufend<(int)bufsize)
    if (buffer[bufend-1]=='\r' || buffer[bufend-1]=='\n')
      break;
  if (bufend == bufpos)
    goteof = true;
  return goteof;
}

size_t 
ParsingByteStream::read(void *buf, size_t size)
{
  if (size < 1)
    return 0;
  if (bufend == bufpos) 
    {
      if (size >= bufsize)
        return bs.read(buf, size);
      if (eof())
        return 0;
    }
  if (bufpos + (int)size > bufend)
    size = bufend - bufpos;
  memcpy(buf, buffer+bufpos, size);
  bufpos += size;
  return size;
}

size_t 
ParsingByteStream::write(const void *, size_t )
{
  G_THROW("Cannot write() into a ParsingByteStream");
  return 0;
}

long int
ParsingByteStream::tell() const
{ 
  G_THROW("Cannot tell() a ParsingByteStream");
  return 0;
}

int 
ParsingByteStream::getbom(int c)
{
  int i = 0;
  while (c == utf8bom[i++])
    {
      if (i >= 3)
        i = 0;
      if (bufpos < bufend || !eof())
        c = buffer[bufpos++];
    }
  while (--i > 0)
    {
      unget(c);
      c = utf8bom[i-1];
    }
  return c;
}

inline int 
ParsingByteStream::get() // like getc() skipping bom.
{
  int c = EOF;
  if (bufpos < bufend || !eof())
    c = buffer[bufpos++];
  if (c == utf8bom[0])
    return getbom(c);
  return c;
}


int  
ParsingByteStream::unget(int c) // like ungetc()
{
  if (bufpos > 0 && c != EOF) 
    return buffer[--bufpos] = (unsigned char)c;
  return EOF;
}

int
ParsingByteStream::get_spaces(bool skipseparator)
{
   int c = get();
   while (c==' ' || c=='\t' || c=='\r' 
          || c=='\n' || c=='#' || c==';' )
     {
       if (c == '#')
         do { c=get(); } while (c!=EOF && c!='\n' && c!='\r');
       if (!skipseparator && (c=='\n' || c=='\r' || c==';'))
         break;
       c = get();
     }
   return c;
}
  
const char *
ParsingByteStream::get_error_context(int c)
{
  static char buffer[22];
  unget(c);
  int len = read((void*)buffer, sizeof(buffer)-1);
  buffer[(len>0)?len:0] = 0;
  for (int i=0; i<len; i++)
    if (buffer[i]=='\n')
      buffer[i] = 0;
  return buffer;
}

GUTF8String
ParsingByteStream::get_token(bool skipseparator, bool compat)
{
  GUTF8String str;
  int c = get_spaces(skipseparator);
  if (c == EOF)
    {
      return str;
    }
  if (!skipseparator && (c=='\n' || c=='\r' || c==';'))
    {
      unget(c);
      return str;
    }
  if (c != '\"' && c != '\'') 
    {
      while (c!=' ' && c!='\t' && c!='\r' && c!=';'
             && c!='\n' && c!='#' && c!=EOF)
        {
          str += c;
          c = get();
        }
      unget(c);
    }
  else 
    {
      int delim = c;
      c = get();
      while (c != delim && c!=EOF) 
        {
          if (c == '\\') 
            {
              c = get();
              if (compat && c!='\"')
                {
                  str += '\\';
                }
              else if (c>='0' && c<='7')
                {
                  int x = 0;
                  { // extra nesting for windows
                    for (int i=0; i<3 && c>='0' && c<='7'; i++) 
                    {
                      x = x * 8 + c - '0';
                      c = get();
                    }
                  }
                  unget(c);
                  c = x;
                }
              else 
                {
                  const char *tr1 = "tnrbfva";
                  const char *tr2 = "\t\n\r\b\f\013\007";
                  { // extra nesting for windows
                    for (int i=0; tr1[i]; i++)
                    {
                      if (c == tr1[i])
                        c = tr2[i];
                    }
                  }
                }
            }
          if (c != EOF)
            str += c;
          c = get();
        }
    }
  return str;
}

