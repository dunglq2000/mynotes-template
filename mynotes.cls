\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mynotes}[2023/11/03 LaTeX Template for Notes]
\LoadClass[oneside]{book}

\RequirePackage[fontsize=14pt]{fontsize}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Language part
\RequirePackage{fontspec}
%\RequirePackage{polyglossia}
%\setdefaultlanguage{vietnamese}
%\setotherlanguages{english,russian}
%\SetLanguageKeys{vietnamese}{indentfirst=true}

\setmainfont{CMU Serif}
\setsansfont{CMU Sans Serif}
%\setmonofont{FreeMono}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layout part
\RequirePackage{geometry}
\geometry{a4paper,top=2cm,left=2cm,bottom=2cm,right=2cm,
        headheight=110pt}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bibliography part
\RequirePackage[
    parentracker=true,
    backend=biber,
    hyperref=auto,
    language=russian,
    autolang=other,
    citestyle=gost-numeric,
    defernumbers=true,
    bibstyle=gost-footnote,
    sorting=none
]{biblatex}

\defbibenvironment{gostbibliography}
  {\list
     {\printtext[labelnumberwidth]{%
	  \printfield{labelprefix}%
	  \printfield{labelnumber}}}
     {\toggletrue{bbx:gostbibliography}%
      \renewcommand*{\revsdnamepunct}{\addcomma}%
      \renewcommand*{\labelnamepunct}{\addperiod\space}%
      %\renewcommand*{\bibname}{\hspace{1.27cm}Список используемых литератур}
      %\setlength{\leftmargin}{\bibhang}%
      %\setlength{\itemindent}{\leftmargin}%
      %\setlength{\itemsep}{\bibitemsep}%
      %\setlength{\parsep}{\bibparsep}
      \setlength{\leftmargin}{0cm}
      \setlength{\itemindent}{1.27cm}
      }
    }
  {\endlist}
  {\item}

%\DeclareFieldFormat{urldate}{(дата обращения #1)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Environment part
\RequirePackage[framemethod=tikz]{mdframed}
\RequirePackage{xcolor}
\RequirePackage{amsthm}
\usetikzlibrary{shadows}


\tikzset{titregris/.style =
    {draw=gray, thick, fill=white, shading = definitiontitle,
    text=gray, rectangle, rounded corners, right,
    minimum height=.7cm}}

% Definition
\pgfdeclarehorizontalshading{definitionbackground}{100bp}
    {color(0bp)=(magenta!40); color(100bp)=(black!5)}
\pgfdeclarehorizontalshading{definitiontitle}{100bp}
    {color(0bp)=(blue!40); color(100bp)=(black!5)}

\newcounter{definition}[chapter]
\newenvironment{definition}[1][]{
    \stepcounter{definition}
    \ifstrempty{#1}
    {
        \mdfsetup{
            frametitle={
                \tikz[baseline=(current bounding box.east),outer sep=0pt]
                \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=definitiontitle]
                {\strut \color{black}{Định nghĩa}~\thedefinition};
            },
        }
    }
    {
        \mdfsetup{
            frametitle={
                \tikz[baseline=(current bounding box.east),outer sep=0pt]
                \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=definitiontitle]
                {\strut \color{black}{Định nghĩa}~\thedefinition.~#1};
            },
        }
    }
    \mdfsetup{innertopmargin=10pt,linecolor=blue!20,%
        linewidth=1pt,topline=true,%
        frametitleaboveskip=\dimexpr-\ht\strutbox\relax,
        roundcorner=10pt,
        nobreak=true,
        apptotikzsetting={%
            \tikzset{mdfbackground/.append style={shading=definitionbackground}}
        }
    }
    \begin{mdframed}[]\relax
}{\end{mdframed}}

% Remark
\pgfdeclarehorizontalshading{remarkbackground}{100bp}
    {color(0bp)=(orange!40); color(100bp)=(black!5)}
\pgfdeclarehorizontalshading{remarktitle}{100bp}
    {color(0bp)=(violet!40); color(100bp)=(black!5)}

\newcounter{remark}[chapter]
\newenvironment{remark}[1][]{
    \stepcounter{remark}
    \mdfsetup{innertopmargin=10pt,linecolor=blue!20,%
        linewidth=1pt,topline=true,%
        frametitle={
            \tikz[baseline=(current bounding box.east),outer sep=0pt]
            \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=remarktitle]
            {\strut \color{black}{Nhận xét}~\theremark};
        },
        frametitleaboveskip=\dimexpr-\ht\strutbox\relax,
        roundcorner=10pt,
        nobreak=true,
        apptotikzsetting={%
            \tikzset{mdfbackground/.append style={shading=remarkbackground}}
        }
    }
    \begin{mdframed}[]\relax
}{\end{mdframed}}

% Theorem
\pgfdeclarehorizontalshading{theorembackground}{100bp}
    {color(0bp)=(yellow!40); color(100bp)=(black!5)}
\pgfdeclarehorizontalshading{theoremtitle}{100bp}
    {color(0bp)=(green!40); color(100bp)=(black!5)}

\newcounter{theorem}[chapter]
\newenvironment{theorem}[1][]{
    \stepcounter{theorem}
    \ifstrempty{#1}
    {
        \mdfsetup{
            frametitle={
                \tikz[baseline=(current bounding box.east),outer sep=0pt]
                \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=theoremtitle]
                {\strut \color{black}{Định lí}~\thetheorem};
            },
        }
    }
    {
        \mdfsetup{
            frametitle={
                \tikz[baseline=(current bounding box.east),outer sep=0pt]
                \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=theoremtitle]
                {\strut \color{black}{Định lí~\thetheorem.~#1}};
            }
        }
    }
    \mdfsetup{innertopmargin=10pt,linecolor=blue!20,%
        linewidth=1pt,topline=true,%
        frametitleaboveskip=\dimexpr-\ht\strutbox\relax,
        roundcorner=10pt,
        nobreak=true,
        apptotikzsetting={%
            \tikzset{mdfbackground/.append style={shading=theorembackground}}
        }
    }
    \begin{mdframed}[]\relax
}{\end{mdframed}}

% Lemma
\pgfdeclarehorizontalshading{lemmabackground}{100bp}
    {color(0bp)=(teal!40); color(100bp)=(black!5)}
\pgfdeclarehorizontalshading{lemmatitle}{100bp}
    {color(0bp)=(blue!40); color(100bp)=(black!5)}

\newcounter{lemma}[chapter]
\newenvironment{lemma}[1][]{
    \stepcounter{lemma}
    \ifstrempty{#1}
    {
        \mdfsetup{
            frametitle={
                \tikz[baseline=(current bounding box.east),outer sep=0pt]
                \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=lemmatitle]
                {\strut \color{black}{Bổ đề}~\thetheorem};
            },
        }
    }
    {
        \mdfsetup{
            frametitle={
                \tikz[baseline=(current bounding box.east),outer sep=0pt]
                \node[line width=1pt,anchor=east,rectangle,rounded corners,draw=blue!20,fill=white,shading=lemmatitle]
                {\strut \color{black}{Bổ đề~\thetheorem.~#1}};
            }
        }
    }
    \mdfsetup{innertopmargin=10pt,linecolor=blue!20,%
            linewidth=1pt,topline=true,%
            frametitleaboveskip=\dimexpr-\ht\strutbox\relax,
            roundcorner=10pt,
            nobreak=true,
            apptotikzsetting={%
            \tikzset{mdfbackground/.append style={shading=lemmabackground}}
        }
    }
    \begin{mdframed}[]\relax
}{\end{mdframed}}

\theoremstyle{definition}
\newtheorem{example}{Ví dụ}
\newtheorem{corollary}{Hệ quả}

\renewcommand*{\proofname}{Chứng minh}
\counterwithin*{example}{chapter}

\renewcommand*{\chaptername}{Chương}
\renewcommand*{\partname}{Phần}
\renewcommand*{\appendixname}{Phụ lục}
\renewcommand*{\contentsname}{Mục lục}
\renewcommand*{\figurename}{Hình}
\renewcommand*{\tablename}{Bảng}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{enumitem}
\RequirePackage{hyperref}
\RequirePackage{graphicx}
\RequirePackage{indentfirst}

\RequirePackage{caption}
\RequirePackage{subcaption}

\DeclareCaptionLabelSeparator{bar}{.\space}
\captionsetup[figure]{labelsep=bar}
\captionsetup[table]{labelsep=bar}
\captionsetup[sub]{font=normalsize}