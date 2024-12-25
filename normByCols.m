function out=normByCols(in)
out=sqrt(sum(in.*conj(in),1));