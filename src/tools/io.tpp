/* template implementation file for io */

template <class T>
T* io::RootFileReader::read(TString name) const
{
   if (intPath_!="") name=intPath_+"/"+name;
   T* h=(T*)file_->Get(name);
   if (!h){
      debug_io<<TString::Format("did not find '%s' in '%s'",(name).Data(),fName_.Data());
      throw;
   }
   return h;
}

