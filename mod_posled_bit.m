clear
clc
%%%%%%%%%%%Моделирование последовательности бит%%%%%%%%%%%
%Номер варианта 3
j=sqrt(-1);
p1=0.35;
Nbit=390*30*5;
for i = 1:Nbit
	b = rand; % генерация числа, равномерно распределённого от 0 до 1
	if(b>p1) 
		msg(i) = 1; % присваивание 1
	else
msg(i) = 0; % присваивание 0
    end
end


%%%%%%%%%%%Моделирование кодера%%%%%%%%%%%
n_b=31;
k_b=26;
n_c=16;
k_c=6;

pol = cyclpoly(n_b,k_b);
parmat = cyclgen(n_b,pol);
genmat = gen2par(parmat);
genpoly = cyclpoly(n_c,k_c);
trellis = poly2trellis(3,[6 7]);

code1=encode(msg,n_b,k_b,'linear/binary',genmat);%Блочный кодер
code2=encode(msg,n_c,k_c,'cyclic/binary',genpoly);%Циклический кодер
code3=convenc(msg,trellis);%Свёрточный кодер
cage=poly2trellis(3,[6 7]);%Решётка для свёрточного кода

%decode1=decode(code1,n_b,k_b,'linear/binary',genmat);%Блочный декодер
%decode2=decode(code2,n_c,k_c,'cyclic/binary',genpoly);%Циклический декодер
%tbdepth=5%A rate 1/2 code has a traceback depth of 5
%decode3=vitdec(code3,trellis,tbdepth,'trunc','hard')%Свёрточный декодер

%%%%%%%%%%%Моделирование перемежителя%%%%%%%%%%%
Nrows=10;
Ncols=3;
i1=1;
i2=Nrows*Ncols;
for i = 1:length(code1)/(Nrows*Ncols)
    op1(i1:i2)=matintrlv(code1(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

i1=1;
i2=Nrows*Ncols;
for i = 1:length(code2)/(Nrows*Ncols)
    op2(i1:i2)=matintrlv(code2(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

i1=1;
i2=Nrows*Ncols;
for i = 1:length(code3)/(Nrows*Ncols)
    op3(i1:i2)=matintrlv(code3(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end


%%%%%%%%%%%Модуляция QAM32%%%%%%%%%%%
M=32;
k=log2(M);%Бит на символ
QAM = qammod(op1.',M,'gray', InputType='bit', UnitAveragePower=true);
scatterplot(QAM)

%%%%%%%%%%%Моделирование канала связи%%%%%%%%%%%
%Гауссовский канал
awgnchan = comm.AWGNChannel('SNR',15,'BitsPerSymbol',k);
QAM_noise = awgnchan(QAM);
scatterplot(QAM_noise)

%%%%%%%%%%%Демодуляция QAM32%%%%%%%%%%%
QAM_demod = qamdemod(QAM.*exp(-1j*pi/M), M,'gray', OutputType='bit',UnitAveragePower=true);
isequal(QAM_demod.',op1)


%%%%%%%%%%%Моделирование деперемежителя%%%%%%%%%%%
Nrows=10;
Ncols=3;
i1=1;
i2=Nrows*Ncols;
for i = 1:length(op1)/(Nrows*Ncols)
    dop1(i1:i2)=matdeintrlv(op1(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

i1=1;
i2=Nrows*Ncols;
for i = 1:length(op2)/(Nrows*Ncols)
    dop2(i1:i2)=matdeintrlv(op2(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

i1=1;
i2=Nrows*Ncols;
for i = 1:length(op3)/(Nrows*Ncols)
    dop3(i1:i2)=matdeintrlv(op3(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

%%%%%%%%%%%Моделирование декодера%%%%%%%%%%%
decode1=decode(dop1,n_b,k_b,'linear/binary',genmat);%Блочный декодер
decode2=decode(dop2,n_c,k_c,'cyclic/binary',genpoly);%Циклический декодер
tbdepth=5;%A rate 1/2 code has a traceback depth of 5
decode3=vitdec(dop3,trellis,tbdepth,'trunc','hard');%Свёрточный декодер

%Проверка на ошибки
if (decode1==msg) & (decode2==msg) & (decode3==msg)
    disp('Ошибок нет')
else
    disp('НООООУ')
end
