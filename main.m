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
QAM_b = qammod(op1.',M,'gray', InputType='bit', UnitAveragePower=true);
QAM_c = qammod(op2.',M,'gray', InputType='bit', UnitAveragePower=true);
QAM_s = qammod(op3.',M,'gray', InputType='bit', UnitAveragePower=true);
%scatterplot(QAM)

%%%%%%%%%%%Моделирование канала связи%%%%%%%%%%%
%Гауссовский канал
for dB = [0 5 10 15 20 25 30]
    disp(dB)
awgn=comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',dB);
%awgnchan = comm.AWGNChannel('SNR',dB,'BitsPerSymbol',k);
QAM_b_noise = awgn(QAM_b);
QAM_c_noise = awgn(QAM_c);
QAM_s_noise = awgn(QAM_s);
%scatterplot(QAM_b_noise)

%%%%%%%%%%%Демодуляция QAM32%%%%%%%%%%%
QAM_b_demod = qamdemod(QAM_b_noise.*exp(-1j*pi/M), M,'gray', OutputType='bit',UnitAveragePower=true);
QAM_c_demod = qamdemod(QAM_c_noise.*exp(-1j*pi/M), M,'gray', OutputType='bit',UnitAveragePower=true);
QAM_s_demod = qamdemod(QAM_s_noise.*exp(-1j*pi/M), M,'gray', OutputType='bit',UnitAveragePower=true);
%isequal(QAM_demod.',op1)


%%%%%%%%%%%Моделирование деперемежителя%%%%%%%%%%%
Nrows=10;
Ncols=3;
i1=1;
i2=Nrows*Ncols;
for i = 1:length(QAM_b_demod)/(Nrows*Ncols)
    dop1(i1:i2)=matdeintrlv(QAM_b_demod(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

i1=1;
i2=Nrows*Ncols;
for i = 1:length(QAM_c_demod)/(Nrows*Ncols)
    dop2(i1:i2)=matdeintrlv(QAM_c_demod(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

i1=1;
i2=Nrows*Ncols;
for i = 1:length(QAM_s_demod)/(Nrows*Ncols)
    dop3(i1:i2)=matdeintrlv(QAM_s_demod(i1:i2),Nrows,Ncols);
    i1=i1+Nrows*Ncols;
    i2=i2+Nrows*Ncols;
end

%%%%%%%%%%%Моделирование декодера%%%%%%%%%%%
decode1=decode(dop1,n_b,k_b,'linear/binary',genmat);%Блочный декодер
decode2=decode(dop2,n_c,k_c,'cyclic/binary',genpoly);%Циклический декодер
tbdepth=5;%A rate 1/2 code has a traceback depth of 5
decode3=vitdec(dop3,trellis,tbdepth,'trunc','hard');%Свёрточный декодер

%isequal(decode1,msg)
%isequal(decode2,msg)
%isequal(decode3,msg)


%Сравнение
er_count_b=0;
for i = 1:length(msg)
    if (decode1(i)~=msg(i))
        er_count_b=er_count_b+1;
    end
end

Po_b(dB+1)=er_count_b/length(msg)*100;

er_count_c=0;
for i = 1:length(msg)
    if (decode2(i)~=msg(i))
        er_count_c=er_count_c+1;
    end
end

Po_c(dB+1)=er_count_c/length(msg)*100;

er_count_s=0;
for i = 1:length(msg)
    if (decode3(i)~=msg(i))
        er_count_s=er_count_s+1;
    end
end

Po_s(dB+1)=er_count_s/length(msg)*100;

end
dB = [0 5 10 15 20 25 30];
plot(dB,Po_b(dB+1),dB,Po_c(dB+1),dB,Po_s(dB+1))
title('Pош от ОСШ. КАМ32. Гауссовский канал')
xlabel('дБ')
ylabel('Pош')
legend('Блочный', 'Циклический', 'Свёрточный')
grid