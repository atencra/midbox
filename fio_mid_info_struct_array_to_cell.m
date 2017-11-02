

indfirst = [1 17 31 50 64 83 92 102 129 140 151 165 183 211 226];

for i = 1:length(indfirst)

if ( i < length(indfirst) )
	fiocell{i} = fiopos( indfirst(i):indfirst(i+1)-1 );
	midcell{i} = midpos( indfirst(i):indfirst(i+1)-1 );
	infocell{i} = infopos( indfirst(i):indfirst(i+1)-1 );

else
	fiocell{i} = fiopos( indfirst(i):end );
	midcell{i} = midpos( indfirst(i):end );
	infocell{i} = infopos( indfirst(i):end );
end

end





