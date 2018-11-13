function []=rd_gen(n,m)

	nn = n+m-1;
	mm = n*m;
	A=zeros(nn,mm);
	for i = 1:m
		for j = 1:n
			A(i,i*n+j)=1;
			A(j+m,i*n+j)=1;
		end
	end

	fa = fopen(['tpA',num2str(nn),'x',num2str(mm),'.mtx'],'w+');
	fprintf(fa,'%%MatrixMarket matrix array real general\r\n');
	fprintf(fa,'%d %d\r\n',nn,mm);
	fprintf(fa,'%f\r\n',A);
	fclose(fa);

	b=zeros(m+n-1,1);
	b(1:m)=ones(m,1)/m;
	b(m+1:n+m-1)=ones(n-1,1)/n;
	fb = fopen(['tpb',num2str(nn),'x',num2str(mm),'.mtx'],'w+');
	fprintf(fb,'%%MatrixMarket matrix array real general\r\n');
	fprintf(fb,'%d %d\r\n',nn,1);
	fprintf(fb,'%f\r\n',b);
	fclose(fb);


	fc = fopen(['tpc',num2str(nn),'x',num2str(mm),'.mtx'],'w+');
	fprintf(fc,'%%MatrixMarket matrix array real general\r\n');
	fprintf(fc,'%d %d\r\n',mm,1);
	fprintf(fc,'%f\r\n',rand(mm,1));
	fclose(fc);

end
