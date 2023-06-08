% codigoFuenteLibroR2
% Copyright (C) 2020  J.M. Perez Zerpa, et. al.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

\def\filename{textbook.cls}
\def\fileversion{v0.1}
\def\filedate{2019/06/28}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{textbook}[\filedate\ \fileversion\ LaTeX Class for writing textbooks in english or spanish.]
%\LoadClass[oneside,openright,10pt]{book}
\LoadClass[oneside,openright,10pt]{book}
%
% --- packages ---
%
%\newcommand\datoscurso[2]{%
%  \global\def\@nombrecurso{#1}%
%  \global\def\@codigocurso{#2}%
%}
%%
%\newcommand\tipodeprueba[1]{%
%  \global\def\@tipodeprueba{#1}%
%}
%%
%\newcommand\universidad[1]{%
%	\global\def\@universidad{#1}%
%}
%%
%\newcommand\carrera[1]{%
%	\global\def\@carrera{#1}%
%}
%

\newcounter{parte}
%
\newcommand{\parte}{
	\stepcounter{parte}
	\vskip 2mm
	\par\noindent{\alph{parte})}\;
}
%
\newcounter{ejercicio}
\setcounter{ejercicio}{0}
\newcommand{\ejercicio}{
	\setcounter{parte}{0}
	\stepcounter{ejercicio}
	\vskip 6mm
	\par\noindent{\large{\bf Ejercicio \arabic{chapter}.\arabic{ejercicio}}}\;
	%    \vskip 4mm
}
%\counterwithin{ejercicio}{chapter}
%
\newcounter{pregunta}
\setcounter{pregunta}{0}
\newcommand{\pregunta}{
	\setcounter{parte}{0}
	\stepcounter{pregunta}
	\vskip 4mm
	\par\noindent{\large{\bf Pregunta \arabic{pregunta}}}\;
}
