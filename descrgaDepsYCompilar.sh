wget https://raw.githubusercontent.com/jorgepz/svg2pdf_tex/master/svg2pdftex.sh
mv svg2pdftex.sh latex/

cd latex
chmod +x svg2pdftex.sh
./svg2pdftex.sh 1
cd ..


wget https://gist.githubusercontent.com/jorgepz/9b3050d12ce9e79cdeac415b3845c715/raw/62f9e19cffde0e1521ba59d0144a9f3038af7cb6/definitions.tex 
mv definitions.tex latex/texs/

#cd latex && latexmk -pdf libroResMat2.tex && cd ..
cd latex && pdflatex libroResMat2.tex && cd ..
