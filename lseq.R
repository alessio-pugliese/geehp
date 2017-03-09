lseq<-function(from,to,length.out){
  exp(seq(log(from), log(to), length.out = length.out))
}