\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{math-mech-sci}[2025/04/11 v0.2.1 Dmitry Luciv]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

\LoadClass[10pt,a5paper]{article}

\RequirePackage{ifxetex}

% AMS packages want to be first...
\RequirePackage{amsmath}

\ifxetex

    \RequirePackage{mathspec}
    \RequirePackage{fontspec}
    \RequirePackage{xltxtra}
    \RequirePackage{xecyr}
    \defaultfontfeatures{Mapping=tex-text,Ligatures={TeX}}

    % Windows fonts
    \iffalse
        \def\MyRomanFont#1{\expandafter#1[]{Times New Roman}}
        \def\MyMonoFont#1{\expandafter#1[]{Courier New}}
        \def\MySansFont#1{\expandafter#1{Arial}}
    \else

        \def\MyRomanFont#1{\expandafter#1{Tempora}[
                Path = fonts/Tempora/,
                Extension=.otf,
                UprightFont=*-Regular,
                ItalicFont=*-Italic,
                BoldFont=*-Bold,
                BoldItalicFont=*-BoldItalic,
                Kerning=On
            ]}

        \def\MySansFont#1{\expandafter#1{LiberationSans}[
                Path = fonts/LiberationSans/,
                Extension=.ttf,
                UprightFont=*-Regular,
                ItalicFont=*-Italic,
                BoldFont=*-Bold,
                BoldItalicFont=*-BoldItalic,
                Kerning=On
            ]}

        \def\MyMonoFont#1{\expandafter#1{IosevkaCustom}[
                Scale = 0.9,
                Path = fonts/Iosevka/,
                Extension=.ttf,
                UprightFont=*-Regular,
                ItalicFont=*-Italic,
                BoldFont=*-Bold,
                BoldItalicFont=*-BoldItalic
            ]}

    \fi

    \MyRomanFont{\setmainfont}
    \MyRomanFont{\setromanfont}
    \MyRomanFont{\newfontfamily{\cyrillicfont}}
    \MyMonoFont{\setmonofont}
    \MyMonoFont{\newfontfamily{\cyrillicfonttt}}
    \MySansFont{\setsansfont}
    \MySansFont{\newfontfamily{\cyrillicfontsf}}

    \RequirePackage{polyglossia}
    \setdefaultlanguage{russian}
    \setotherlanguages{english}

\else

    \RequirePackage{mmap}
    \RequirePackage[T2A]{fontenc}
    \RequirePackage[utf8]{inputenc}
    \RequirePackage[english,russian]{babel}
    \selectlanguage{russian}

\fi

\RequirePackage[unicode,pdfa]{hyperref}
\RequirePackage{booktabs}
\RequirePackage{indentfirst}
\RequirePackage{titlesec}
\RequirePackage{graphicx}
\RequirePackage[table,xcdraw]{xcolor}
\RequirePackage{listings}
%\RequirePackage{minted}
\RequirePackage{makecell}
\RequirePackage{amscd}
\RequirePackage{amssymb}
\RequirePackage{accents}
\RequirePackage[author={Мат-мех. Наука}]{pdfcomment}

\RequirePackage[top=23mm,left=17mm,right=17mm,bottom=17mm]{geometry}

\RequirePackage[font=small,skip=0pt]{caption}

\pagestyle{empty}

\renewenvironment{abstract}
{\normalsize
    \begin{center}
        \bfseries\large \abstractname \normalsize\vspace{-.5em}
    \end{center}
    \list{}{
        \setlength{\leftmargin}{10mm}%
        \setlength{\rightmargin}{\leftmargin}%
        \setlength{\listparindent}{\parindent}%
    }%
    \item\relax}
{\endlist}

% remove numbering
\setcounter{secnumdepth}{-1}

\titleformat{\section}
{\center\normalfont\fontsize{12}{17}\bfseries}
{\thesection}
{2pt}
{}

\titleformat{\subsection}
{\center\normalfont\fontsize{11}{17}\bfseries\itshape}
{\thesubsection}
{2pt}
{}

\renewcommand{\thanks}[1]{\normalfont\mdseries\footnote{\footnotesize\textup{#1}}}

\lstset{ %
    basicstyle=\normalsize\ttfamily,  % the size of the fonts that are used for the code
    showstringspaces=false,         % underline spaces within strings
    tabsize=2,                      % sets default tabsize to 2 spaces
    breaklines=true,                % sets automatic line breaking
    breakatwhitespace=false,        % sets if automatic breaks should only happen at whitespace
    escapeinside={\%*}{*)},          % if you want to add a comment within your code
    extendedchars=\true,
    inputencoding=utf8,
    keywordstyle=\color{blue}\bfseries\ttfamily,
    stringstyle=\color{red},
    commentstyle=\color{brown},
    xleftmargin=5mm
}

\newcommand{\thd}[1]{\thead[c]{\normalsize\bfseries #1}}

% http://tex.stackexchange.com/a/232874
\setlength{\tabcolsep}{0.5em} % for the horizontal padding
{\renewcommand{\arraystretch}{1.2}% for the vertical padding

\sloppy

% authors metadata

\RequirePackage{hyperxmp}
\RequirePackage{etoolbox}

\newbool{firstauthor}
\booltrue{firstauthor}
\newcommand{\defauthor}[4][]{%
    \ifbool{firstauthor}{%
        \appto\atext{#2, #3 \href{mailto:#4}{\texttt{#4}}#1}%
        \appto\ameta{#2}%
        \appto\xmeta{\xmpquote{#2}}%
        \boolfalse{firstauthor}
    }{%
        \appto\atext{,\\ #2, #3 \href{mailto:#4}{\texttt{#4}}#1}%
        \appto\ameta{ \and #2}%
        \appto\xmeta{ \and \xmpquote{#2}}%
    }%
}

\def\titlethanks{}
\newcommand{\deftitle}[2][]{%
    \def\titlethanks{#1}
    \title{#2}
}

\renewcommand\maketitle{%
    \author{\ameta}%
    \patchcmd{\xmeta}{, }{\xmpcomma\ }{\author{\xmeta}}{}%
    \hypersetup{%
        pdfauthor={\xmeta}%
    }%
    \begin{center}%
        \vspace*{-12mm}
        {\Large\bfseries \@title\titlethanks \par\vspace{3mm}}%
        {\small\atext}\ifxetex\else
            \pdfcomment{Пожалуйста используйте XeLaTeX}
        \fi
    \end{center}\par%
}

\iffalse


    \patchcmd{\xmeta}{, }{\xmpcomma\ }{\author{\xmeta}}{}%

    \author{\xmeta}
    \hypersetup{%
        pdfauthor={\xmeta}%
    }%
    \show\xmeta%
\fi

\endinput
