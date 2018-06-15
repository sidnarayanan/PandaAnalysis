#!/usr/bin/env python

from sys import argv,exit
from os import system, path
from re import sub, match
import argparse
import copy

parser = argparse.ArgumentParser(description='build object from configuration')
parser.add_argument('--cfg',type=str)
args = parser.parse_args()


suffixes = { 'float':'F', 
             'int':'I',
             'ULong64_t':'l',
           }


def get_word(line, predef):
    line = sub('.*=','',line)
    if ':' in line:
        return predef[line]
    else:
        return line

class Word(object):
    __slots__ = ['name', 'expr']
    def __init__(self, line, predef):
        # word must look like name expr 
        ll = line.split()
        assert len(ll) == 2, "Bad definition: "+line
        self.name = ll[0]
        self.expr = ll[1]
        for n,x in predef.iteritems():
            if n in self.expr:
                self.expr = self.expr.replace(n, str(x))
    def __str__(self):
        return str(self.expr.replace('const:',''))
    def to_tstring(self):
        if self.expr.isdigit():
            return 'TString("%s")'%self.expr 
        else:  # assume it's a bare string or C++ expr
            return 'TString(%s)'%self.expr 
    def write_header(self):
        r = []
        if self.name.startswith('const:'):
            r.append('#define %s %s'%(self.name.replace('const:',''),self.expr))
        return r
    def write_def(self):
        return []

class Shift(Word):
    def __init__(self, line, predef):
        super(Shift, self).__init__(line, predef)
        tmp_expr = self.expr.replace('[','').replace(']','')
        self.expr = []
        self.expr += [('' if not x.strip() else '_'+x.strip()) for x in tmp_expr.split(',')]
    def __str__(self):
        return '{' + ','.join(['"%s"'%x for x in self.expr]) + '}'
    @property
    def cname(self):
        return self.name.replace('shift:','')
    def is_subset(self):
        # this is a subset of another definition  
        return bool(match('.*\[[0-9]*\]$', self.name))
    def write_header(self):
        if self.is_subset():
            return []
        enum = []
        enum.append('enum class shift%s { '%self.cname)
        for e in self.expr:
            if e == '':
                enum.append('  kNominal,')
            else:
                enum.append('  k%s,'%(e[1:]))
        enum.append('  N')
        enum.append('}; ')
        enum.append('TString %sName(shift%s);'%(self.cname,self.cname))
        return enum
    def write_def(self):
        if self.is_subset():
            return []
        fn =     ['TString %sName(shift%s e) {'%(self.cname, self.cname)]
        fn.append('  switch(e) {')
        for e in self.expr:
            if e == '':
                e = '_Nominal'
            e = e.replace('_','')
            fn.append('    case shift%s::k%s : return "%s";'%(self.cname, e, e))
        fn.append('    default : return "N";') 
        fn.append('  }')
        fn.append('}')
        return fn


class Block(object):
    __slots__ = ['enclosing', 'words']
    def __init__(self, line, predef, enclosing):
        self.enclosing = enclosing
        self.words = {}
        if line == '':
            return
        if enclosing:
            self.words = copy.deepcopy(enclosing.words)
        line = line.replace('block:','')
        for w in [x.strip() for x in line.split(',')]:
            word = get_word(w, predef)
            n = w.split('=')[0].strip()
            if n not in self.words:
                self.words[n] = []
            self.words[n].append(word)

class Branch(object):
    __slots__ = ['name', 'array', 'ctype', 'N', 'filled', 'tree_size', 'shift', 'default']
    def __init__(self, line, block, predef):
        args = line.split()
        self.name = args[0]
        self.filled = block.words.get('filled',[])[:]
        self.shift = block.words.get('shifts',[None])[-1]
        self.tree_size = block.words.get('tree_size',[None])[0]
        if '[' in args[1]:
            self.array = True
            self.ctype = sub('\[[:A-z0-9]*\]','',args[1])
            self.N = sub('.*\[([:A-z0-9]*)\]',r'\1',args[1])
            if self.N.startswith('const:'):
                self.N = predef[self.N]
            if not self.tree_size:
                self.tree_size = self.N
            for a in args[2:]:
                if 'tree_size' in a:
                    self.tree_size = Word(a, predef)
            if type(self.tree_size) == str:
                self.tree_size = Word("tree_size "+self.tree_size, predef)
        else:
            self.array = False
            self.ctype = args[1]
        self._set_default() # default defaults
        for a in args[1:]:
            if a.startswith('filled='):
                self.filled.append(get_word(a, predef))
            elif a.startswith('shifts='):
                self.shift = get_word(a, predef)
            elif a.startswith('default='):
                self.default = get_word(a, predef)
    def _set_default(self):
        if self.name.startswith('sf_'):
            self.default = 1
        elif self.ctype == 'int':
            self.default = 0 
        else:
            self.default = -99;
    def write_public(self):
        s = self.ctype + ' ' + self.name 
        if self.shift:
            s += '[' + str(len(self.shift.expr)) + ']'
        if self.array:
            s += '[' + str(self.N) + ']'
        s += ';'
        return s 
    def write_constructor(self):
        return '' 
    def write_destructor(self):
        return '' 
    def write_reset(self):
        s = self.name 
        if self.shift:
            s += '[iS]'
        if self.array:
            s += '[iA]'
        s += ' = %s;'%(str(self.default))
        return s 
    def write_write(self):
        final_s = []
        if self.shift:
            if self.array:
                for i,e in enumerate(self.shift.expr):
                    s = 'Book("{n}{v}",{n}[{i}],"{n}{v}["+{s}+"]/{t}");'.format(
                            n=self.name,
                            v=e,
                            i=i,
                            s=self.tree_size.to_tstring(),
                            t=suffixes[self.ctype]
                        )
                    final_s.append(s)
            else:
                for i,e in enumerate(self.shift.expr):
                    s = 'Book("{n}{v}",&({n}[{i}]),"{n}{v}/{t}");'.format(
                            n=self.name,
                            v=e,
                            i=i,
                            t=suffixes[self.ctype]
                        )
                    final_s.append(s)
        else:
            if self.array:
                s = 'Book("{n}",{n},"{n}["+{s}+"]/{t}");'.format(
                        n=self.name,
                        s=self.tree_size.to_tstring(),
                        t=suffixes[self.ctype]
                    )
                final_s.append(s)
            else:
                s = 'Book("{n}",&{n},"{n}/{t}");'.format(
                        n=self.name,
                        t=suffixes[self.ctype]
                    )
                final_s.append(s)
        return final_s

def extract_custom(lines):
    custom = {}
    current = None
    for line in lines:
        if 'STARTCUSTOM' in line:
            current = line.split()[-1]
            custom[current] = []
        if current:
            custom[current].append(line.replace('\n',''))
        if 'ENDCUSTOM' in line:
            current = None
    return custom

def get_custom(d, n):
    l = d.get(n, None)
    if l is None:
        l = ['// STARTCUSTOM %s'%n, '// ENDCUSTOM']
    return l 

if __name__ == '__main__':
    cfg_path = args.cfg
    header_path = cfg_path.replace('config','interface').replace('.cfg','.h')
    src_path = cfg_path.replace('config','src').replace('.cfg','.cc')
    if not path.isfile(header_path):
        tmpl_path = '/'.join(header_path.split('/')[:-2])
        tmpl_path += ('bin/treeClassTmpl.h' if len(tmpl_path) == 0 else '/bin/treeClassTmpl.h')
        system('cp -v %s %s'%(tmpl_path, header_path))
    if not path.isfile(src_path):
        tmpl_path = '/'.join(src_path.split('/')[:-2])
        tmpl_path += ('bin/treeClassTmpl.cc' if len(tmpl_path) == 0 else '/bin/treeClassTmpl.cc')
        system('cp -v %s %s'%(tmpl_path, src_path))
    class_name = cfg_path.split('/')[-1].replace('.cfg','')

    for p in [header_path, src_path]:
        system('cp -v {0} {0}.bkp'.format(p))

    # first let's parse the cfg file 
    fcfg = open(cfg_path)
    defined = {}
    branches = []
    default_block = Block('', defined, None)
    current_block = default_block 
    for line in fcfg:
        if current_block is None:
            print 'Malformed block structure!'
            exit(1)
        line = line.strip()
        line = sub('#.*','',line)
        if line == '':
            continue 
        if any([line.startswith(x) for x in ['const:','cond:']]):
            word = Word(line, defined)
            defined[word.name] = word 
        elif line.startswith('shift:'):
            word = Shift(line, defined)
            defined[word.name] = word
        elif line.startswith('block:'):
            current_block = Block(line, defined, current_block)
        elif line.startswith('endblock'):
            current_block = current_block.enclosing
        else:
            b = Branch(line, current_block, defined)
            branches.append(b)

    # now let's write the header
    header_tmpl = list(open(header_path,'r'))
    custom = extract_custom(header_tmpl)

    header_out = ['// THIS FILE IS AUTOGENERATED //']
    header_out += ['#ifndef %s_H'%class_name, '#define %s_H'%class_name] 
    header_out.extend(get_custom(custom, 'HEADER'))
    for _,d in defined.iteritems():
        header_out.extend(d.write_header())

    header_out.append(('\n'.join([
            'class {n} : public genericTree {{',
            '  public:',
            '    {n}();',
            '    ~{n}();',
            '    void WriteTree(TTree* t);',
            '    void Fill() {{ treePtr->Fill(); }}',
            '    void Reset();'
            '    void SetAuxTree(TTree*);'
        ])).format(n=class_name))
    header_out.extend(get_custom(custom, 'PUBLIC'))
    header_out.append('  private:')
    header_out.extend(get_custom(custom, 'PRIVATE'))
    header_out.append('  public:')
    for b in branches:
        header_out.append('  ' + b.write_public())
    header_out.append('};')
    header_out.append('#endif')
    with open(header_path,'w') as fout:
        fout.write('\n'.join(header_out))


    # now we do the src file, same deal as before 
    by_array = {}
    by_shift = {}
    by_shift_array = {}
    by_filled = {}
    for b in branches:
        # group based on shared properties
        if b.filled:
            t = tuple(b.filled)
            if t not in by_filled:
                by_filled[t] = []
            by_filled[t].append(b)
        if b.array and b.shift:
            t = (b.shift, b.N)
            if t not in by_shift_array:
                by_shift_array[t] = []
            by_shift_array[t].append(b)
        elif b.array:
            if b.N not in by_array:
                by_array[b.N] = []
            by_array[b.N].append(b)
        elif b.shift:
            if b.shift not in by_shift:
                by_shift[b.shift] = []
            by_shift[b.shift].append(b)
    src_tmpl = list(open(src_path,'r'))
    custom = extract_custom(src_tmpl)

    src_out = ['// THIS FILE IS AUTOGENERATED //']
    src_out += ['#include "../interface/%s.h"'%class_name] 
    src_out.extend(get_custom(custom, 'INCLUDE'))

    for _,d in defined.iteritems():
        src_out.extend(d.write_def())

    src_out.append('{n}::{n}() {{'.format(n=class_name))
    src_out.extend(get_custom(custom, 'CONSTRUCTOR'))
    src_out.append('  Reset();')
    src_out.append('}')

    src_out.append('{n}::~{n}() {{'.format(n=class_name))
    src_out.extend(get_custom(custom, 'DESTRUCTOR'))
    src_out.append('}')

    src_out.append('void {n}::SetAuxTree(TTree *t) {{'.format(n=class_name))
    src_out.extend(get_custom(custom, 'AUX'))
    src_out.append('}')

    src_out.append('void {n}::Reset() {{'.format(n=class_name))
    src_out.extend(get_custom(custom, 'RESET'))
    for b in branches:
        if not b.array and not b.shift:
            src_out.append('  ' + b.write_reset())
    for N,bs in by_array.iteritems():
        src_out.append('  for (int iA=0; iA!=%s; ++iA) {'%N)
        for b in bs:
            src_out.append('    '+b.write_reset())
        src_out.append('  }')
    for s,bs in by_shift.iteritems():
        N = len(s.expr)
        src_out.append('  for (int iS=0; iS!=%s; ++iS) {'%N)
        for b in bs:
            src_out.append('    '+b.write_reset())
        src_out.append('  }')
    for (s,NA),bs in by_shift_array.iteritems():
        NS = len(s.expr)
        src_out.append('  for (int iS=0; iS!=%s; ++iS) {'%NS)
        src_out.append('    for (int iA=0; iA!=%s; ++iA) {'%NA)
        for b in bs:
            src_out.append('      '+b.write_reset())
        src_out.append('    }')
        src_out.append('  }')
    src_out.append('}')

    src_out.append('void {n}::WriteTree(TTree *t) {{'.format(n=class_name))
    src_out.append('  treePtr = t;')
    src_out.extend(get_custom(custom, 'WRITE'))
    for b in branches:
        if not b.filled:
            src_out.extend(['  ' + x for x in b.write_write()])
    for cond,bs in by_filled.iteritems():
        for i,c in enumerate(cond):
            src_out.append('  ' + ('  ' * i) + 'if (%s) {'%str(c))
        for b in bs:
            src_out.extend(['    ' + ('  ' * i) + x for x in  b.write_write()])
        for i,c in enumerate(cond):
            src_out.append('  ' + ('  ' * (len(cond) - i - 1)) + '}')
    src_out.append('}')
    with open(src_path,'w') as fout:
        fout.write('\n'.join(src_out))



