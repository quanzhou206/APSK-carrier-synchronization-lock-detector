clc;clear;

M = [4 12];
radii = [sqrt(2) 3];

k = log2(16);
numBits = 1024*4; % 传输总bit
N_sym = numBits/k;
epoch_all = 2000;
fs = 16e3; % 采样
Rs = 1e3; % 符号率
mode = 1;% 1 lock;2 unlock
if mode==1
    fd = 0;
else
    fd = rand()*1e3; % 载波频差
end
sps = fs/Rs; % 每符号上采样数
filtlen = 5; % 滤波窗口长度
rolloff = 0.25; % 滚将系数

rrcFilter = rcosdesign(rolloff,filtlen,sps);

dataIn = randi([0 1],numBits,1); % bit数据
dataInMatrix = reshape(dataIn,length(dataIn)/k,k); % 数据转矩阵
dataSymbolsIn = bi2de(dataInMatrix); % 矩阵转16进制
dataMod = apskmod(dataSymbolsIn,M,radii); % qam调制
txFiltSignal = upfirdn(dataMod,rrcFilter,sps,1); % 成型滤波

% t = 0:1/fs:(length(txFiltSignal)-1)/fs;
t = 0:1/fs:(length(dataMod)-1)/fs;
carrier = exp(1j*2*pi*fd*t);
txSignal = dataMod.*carrier.';
% txSignal = carrier.';

EbNo = -5:2:25; % QAM平均Eb/n0 dB

% snr = EbNo + 10*log10(k) - 10*log10(sps);
snr = EbNo + 10*log10(k);

Seg_th = (sqrt(2)+3)/2;
weight_matrix = [2/11,9/11];% 2个圈权重

for epoch = 1:epoch_all
for index = 1:length(snr)
    rxSignal = awgn(txSignal,snr(index),'measured');

%     rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,sps); % 匹配滤波
%     rxFiltSignal = rxFiltSignal(filtlen + 1:end - filtlen); % 消除滤波延时
    rxFiltSignal = rxSignal;
    xMN_w = zeros(1,2);
    xMN_w_num = zeros(1,2);
    xMN = 0;
    for i = 1:length(rxFiltSignal)
        AmFiltSignal = abs(rxFiltSignal(i)); % 获取信号幅度
%         dataclass1 = get_class(abs(dataMod(i)),Seg_th);
        dataclass = get_class(AmFiltSignal,Seg_th); %获取每个数据类别

        if dataclass == 1
            xMN_w(dataclass) = xMN_w(dataclass) + ...
                            real(rxFiltSignal(i).^4)./abs(rxFiltSignal(i)).^4;
        else
            xMN_w(dataclass) = xMN_w(dataclass) + ...
                            real(rxFiltSignal(i).^12)./abs(rxFiltSignal(i))^12;
        end
        xMN = xMN + (real(rxFiltSignal(i)^12)./abs(rxFiltSignal(i))^12);
        xMN_w_num(dataclass) = xMN_w_num(dataclass) + 1;
    end
    dataDetect_c1(epoch,index) = abs(xMN)/N_sym;
    dataDetect_c2(epoch,index) = (abs(xMN_w(2))./xMN_w_num(2));
    dataDetect_w(epoch,index) = (abs(xMN_w)./xMN_w_num)*weight_matrix';
end
end
%% 仿真E sig

for i=1:length(snr)
    ELockDetect_w(1,i) = mean(dataDetect_w(:,i));
    SigLockDetect_w(1,i) = std(dataDetect_w(:,i));
    ELockDetect_c1(1,i) = mean(dataDetect_c1(:,i));
    SigLockDetect_c1(1,i) = std(dataDetect_c1(:,i));
    ELockDetect_c2(1,i) = mean(dataDetect_c2(:,i));
    SigLockDetect_c2(1,i) = std(dataDetect_c2(:,i));
end
% 理论E
SigAm = ([2,9]);
averAm = sum(abs(dataMod).^2)/numel(dataMod);
y = zeros(2,length(EbNo));
y_2 = zeros(2,length(EbNo));
for i=1:2
    EsN0_class = snr +  10*log10(SigAm(i)/averAm);
    ESN0 = 10.^(EsN0_class./10);
    if i==1
        y(i,:) = exp(-4^2./(4*ESN0));
        y_2(i,:) = ((0.5 + exp(-4^2./ESN0)/2) - y(i,:).^2)/(N_sym/4);
    else
        y(i,:) = exp(-12^2./(4*ESN0));
        y_2(i,:) = ((0.5 + exp(-12^2./ESN0)/2) - y(i,:).^2)/(3*N_sym/4);
    end
end
y_aver_MN = exp(-12^2./(4*(10.^(snr./10))));
Sig_aver_MN = (0.5 + exp(-12^2./(10.^(snr./10)))/2) - y_aver_MN.^2;
Sig_aver_MN = sqrt(Sig_aver_MN/N_sym);

Sig_c2 = sqrt(y_2(2,:));

y_aver_w = (weight_matrix)*y;
y_sig_w = sqrt((weight_matrix.^2*y_2));

%
figure;hold on;
plot(snr,ELockDetect_w,'-s','linewidth',1.0);grid on;
plot(snr,y_aver_w,'^','linewidth',1.0);
plot(snr,ELockDetect_c1,'-s','linewidth',1.0);grid on;
% plot(snr,y_aver_MN,'^','linewidth',1.0);
plot(snr,ELockDetect_c2,'-s','linewidth',1.0);grid on;
% plot(snr,y(2,:),'*','linewidth',1.0);
% plot(EbNo + 10*log10(k),ELockDetect,'-*','linewidth',1.0);
xlabel('Es/N0');
ylabel('Detection probability');
legend('simulation','theoretical','Linn','最外圆周');
axis tight

figure;hold on;
plot(snr,SigLockDetect_w,'-^','linewidth',1.0);grid on;
plot(snr,y_sig_w,'-s','linewidth',1.0);grid on;
plot(snr,SigLockDetect_c1,'-^','linewidth',1.0);grid on;
% plot(snr,Sig_aver_MN,'-s','linewidth',1.0);grid on;
plot(snr,SigLockDetect_c2,'-^','linewidth',1.0);grid on;
% plot(snr,Sig_c2,'-s','linewidth',1.0);grid on;
xlabel('Es/N0');
ylabel('Detection sigma');
legend('simulation','theoretical','Linn','最外圆周');
axis tight


%%
function dataclass = get_class(data,Seg_th)
    dataclass = zeros(1,length(data));
    for i=1:length(data)
        Am = abs(data(i));
        if Am <= Seg_th
            dataclass(i) = 1;
        elseif Am > Seg_th
            dataclass(i) = 2;
        end
    end
end