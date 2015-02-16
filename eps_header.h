#ifndef EPS_HEADER
#define EPS_HEADER

#define EPS_HEAD "\
%!PS-Adobe-3.0 EPSF-3.0\n\
%%BoundingBox: 0 0 352 420\n\
/Helvetica findfont\n\
20 scalefont\n\
0 setgray\n\
setfont\n\
\n\
/cL { % w x y\n\
 newpath\n\
 0 setgray\n\
 moveto\n\
 0 rlineto\n\
 1 setlinewidth\n\
 stroke\n\
} def\n\
\n\
/conn { % x0 y0 x1 y1 cx cy\n\
 newpath\n\
 moveto\n\
 lineto\n\
 0 setgray\n\
 stroke\n\
} def\n\
\n\
/dY 16 def\n\
\n\
/exon{ % r, g, b, w, x, y\n\
 newpath\n\
 moveto % x, y\n\
 0 dY rlineto\n\
 0 rlineto % w\n\
 0 dY neg rlineto\n\
 closepath\n\
 gsave\n\
 setrgbcolor fill\n\
 % r g b\n\
 grestore \n\
 0 setgray\n\
 stroke \n\
} def\n\
"


#endif
