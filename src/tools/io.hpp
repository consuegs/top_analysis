#ifndef IO_HPP__
#define IO_HPP__

#include <TFile.h>
#include <TObject.h>
#include <TString.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <filesystem>


#include "tools/gfx.hpp"

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define debug_io io::log*"DEBUG:"*__FILENAME__*":"*__FUNCTION__*":"*__LINE__*": "

namespace io
{
   std::string shellOutput(const char* cmd);

   void ensurePathForFile(TString fileName);
   
   bool fileExists(const std::string&);
   
   std::vector<TString> listAllObjectPaths(TFile* tf, const TString &internal_path, std::vector<TString> &output, bool lastFolder=false);

   /* class for saving Objects (e.g. Canvases) to a root file */
   class RootFileSaver
   {
   public:
      RootFileSaver(TString rootFileName,TString internalPath,bool lumiText=true,bool standardOutputDirectory=true,bool recreate=false);
      ~RootFileSaver();
      void save(TObject const &obj, TString name, bool decorate=true,bool simulation=true,bool addPDF=false) const;
      void save(gfx::SplitCan &obj, TString name,bool simulation=true,bool addPDF=false) const;
      void closeFile() const;
      TString getInternalPath() const;
      TString getFilePath() const;
   private:
      TFile *file_;
      TString fName_,fPath_,intPath_;
      bool bSimulation_,bLumiText_;
   };

   /* class for reading Objects (e.g. histograms) from a root file */
   class RootFileReader
   {
   public:
      RootFileReader(TString rootFileName,TString internalPath="",bool standardOutputDirectory=true,bool useDataBasePath=false);
      ~RootFileReader();
      template <class T> T* read(TString name) const;
      std::vector<TString> listPaths(bool lastFolder=false) const;
      void closeFile() const;
   private:
      TFile *file_;
      TString fName_,fPath_,intPath_;
   };

   /* Logger class. Can be used to print to stdout or to a file.
    * usage:
    *   io::log*"Start"/"Space"*"Nospace"<<"Endnopace";
    *   io::log*"Start"/"Space"*"Nospace">>"Endspace";
    * debug output with filename and linenumber:
    *   debug<<"I am a debug output!";
    */
   class Logger
   {
   public:
      Logger();
      Logger(std::string fileName);
      ~Logger();

      void putVersion();

      template<typename T>
      inline Logger& operator<<(const T& val){
         *ostream_ << val << std::endl;
         return *this;
      }
      template<typename T>
      inline Logger& operator>>(const T& val){
         *ostream_ << " " << val << std::endl;
         return *this;
      }
      template<typename T>
      inline Logger& operator*(const T& val){
         *ostream_ << val;
         return *this;
      }
      template<typename T>
      inline Logger& operator/(const T& val){
         *ostream_ << " " << val;
         return *this;
      }
      // could be used, same precedence as * /
      // template<typename T>
      // inline Logger& operator%(const T& val){
      //    std::cout << val << std::endl;
      //    return *this;
      // }
      void flush() {
         *ostream_<<std::flush;
      }
   private:
      TString fName_;
      std::ostream *ostream_;
   };

   extern Logger log;

} // namespace io


#include "io.tpp" // template implemenatations
#endif /* IO_HPP__ */
