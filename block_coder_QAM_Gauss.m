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

pol = cyclpoly(n_b,k_b);
parmat = cyclgen(n_b,pol);
genmat = gen2par(parmat);

code1=encode(msg,n_b,k_b,'linear/binary',genmat);%Блочный кодер

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

%%%%%%%%%%%Модуляция QAM32%%%%%%%%%%%
M=32;
k=log2(M);%Бит на символ
QAM = qammod(op1.',M,'gray', InputType='bit', UnitAveragePower=true);
%scatterplot(QAM)

%%%%%%%%%%%Моделирование канала связи%%%%%%%%%%%
%Гауссовский канал
dB=30;
awgn=comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',dB);
%awgnchan = comm.AWGNChannel('SNR',dB,'BitsPerSymbol',k);
QAM_noise = awgn(QAM);
%scatterplot(QAM_noise)

%%%%%%%%%%%Демодуляция QAM32%%%%%%%%%%%
QAM_demod = qamdemod(QAM_noise.*exp(-1j*pi/M), M,'gray', OutputType='bit',UnitAveragePower=true);
%isequal(QAM_demod.',op1)


%%%%%%%%%%%Моделирование деперемежителя%%%%%%%%%%%
Nrows=10;
Ncols=3;
i1=1;
i2=Nrows*Ncols;
for i = 1:length(QAM_demod)/(Nrows*Ncols)
    dop1(i1:i2)=matdeintrlv(QAM_demod(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

%%%%%%%%%%%Моделирование декодера%%%%%%%%%%%
decode1=decode(dop1,n_b,k_b,'linear/binary',genmat);%Блочный декодер
isequal(decode1,msg)

%Сравнение
er_count=0;
Po=0;
for i = 1:length(msg)
    if (decode1(i)~=msg(i))
        er_count=er_count+1;
    end
end
er_count
Po=er_count/length(msg)*100