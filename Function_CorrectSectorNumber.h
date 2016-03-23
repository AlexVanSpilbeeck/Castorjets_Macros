int CorrectSectorNumber( int sector ){

  if( sector > 8 ){ return sector - 8; }
  else{ return sector + 8; }

}
