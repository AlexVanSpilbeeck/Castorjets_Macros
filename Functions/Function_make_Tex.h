//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

#include <iostream>
#include <fstream>
using namespace std;

void Four_plots(TString plot0, TString plot1, TString plot2, TString plot3, TString title, TString file){


  ofstream Tex_file;

  Tex_file.open( file , ios::out | ios::app );

  title.ReplaceAll("#eta","$\\eta$");
  title.ReplaceAll("E_{I}", "$E_I$");

  Tex_file 	<< "\\begin{frame}" << endl;
  Tex_file	<< " \\frametitle{" << title << "}" << endl;
  Tex_file	<< "  \\begin{center}" << endl;
  Tex_file	<< "   \\begin{tikzpicture}" << endl;
  Tex_file	<< "    \\node at (0,0) [left] 	{ \\includegraphics[width=.24\\linewidth]{Figures/" 		<< plot0 	<< "} };"	<< endl;
  Tex_file      << "    \\node at (.25\\linewidth,0) [left] { \\includegraphics[width=.24\\linewidth]{Figures/" << plot1 	<< "} };" 	<< endl;
  Tex_file      << "    \\node at (.50\\linewidth,0) [left] { \\includegraphics[width=.24\\linewidth]{Figures/" << plot2 	<< "} };" 	<< endl;
  Tex_file      << "    \\node at (.75\\linewidth,0) [left] { \\includegraphics[width=.24\\linewidth]{Figures/" << plot3 	<< "} };" 	<< endl;
  Tex_file      << "   \\end{tikzpicture}" << endl;
  Tex_file      << "  \\end{center}" << endl;
  Tex_file 	<< "\\end{frame}" << endl;

  Tex_file.close();   
}




void One_plots(TString plot0, TString title, TString file){

  ofstream Tex_file;
  Tex_file.open( file, ios::out | ios::app );

  Tex_file      << "\\begin{frame}" << endl;
  Tex_file      << " \\frametitle{" << title << "}" << endl;
  Tex_file      << "  \\begin{center}" << endl;
  Tex_file      << "   \\begin{tikzpicture}" << endl;
  Tex_file      << "    \\node at (0,0) { \\includegraphics[width=.7\\linewidth]{Figures/" << plot0 << "} };" << endl;
  Tex_file      << "   \\end{tikzpicture}" << endl;
  Tex_file      << "  \\end{center}" << endl;
  Tex_file      << "\\end{frame}" << endl;

  Tex_file.close();
}
