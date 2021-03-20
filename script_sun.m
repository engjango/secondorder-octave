% ==============================================================================
%
% script_sun.m (C) 2021 J. G. Silva <janderson.gomes@ufrr.br>
%
% Este script escolhe aleatoriamente dois arquivos Script_2Ord*.m presentes
% no mesmo diretório e realiza o cálculo, para cada um, dos seguintes 
% parâmetros dos sistemas de segunda ordem, respectivamente:
%
% - Frequência natural (wn); 
% - Fator de amortecimento (zeta); 
% - Função de Transferência [G(s)];
% - Instante de pico (Tp);
% - Ultrapassagem de percentual (UP);
% - Tempo de Acomodação (Ts);
% - Tempo de subida (Tr);
%
% Última modificação: 26/02/2021
%
% ==============================================================================

% PARÂMETROS DE ESCOPO ---------------------------------------------------------

clear all % Limpar todas as variáveis locais e globais definidas pelo usuário e 
          % todas as funções da tabela de símbolos.
close all % Fechar janela(s) de figuras, inclusive todas as figuras com abas 
          % visíveis (HandleVisibility = "on") são fechadas.
clc       % Limpe a tela do terminal e mova o cursor para o canto superior 
          % esquerdo.

pkg load control  % Carrega o pacote "Control" que inclui ferramentas de design 
                  % de sistema de controle auxiliado por computador (CACSD) 
                  % para GNU/Octave, com base na biblioteca SLICOT comprovada.
                  
% CONSTANTES ÚTEIS -------------------------------------------------------------

TOLERANCIA_TEMPO_DE_ACOMODACAO = 2; % Variação na curva (Y) de 10%, 5%, 2% ou 1%

% VARIÁVEIS GLOBAIS ÚTEIS ------------------------------------------------------

fileNames = {}; % Lista com os nomes de todos os arquivos .mat disponíveis

% Sub-struct para definir os parâmetros de um sistema
sistema = struct (  
  "Tp", {-1},                       % Instante de Pico
  "Ts", {-1},                       % Tempo de acomodação
  "UP", {-1},                       % Percentual de ultrapassagem
  "zeta", {-1},                     % Taxa de amortecimento
  "Tr", {-1},                       % Tempo de subida
  "Wn", {-1}                        % Frequência natural
);

% Parâmetros do primeiro sistema
sistemaA = struct (
  "nome_do_arquivo", {"empty"},     % Nome do arquivo carregado na memória
  "vetor_t", {[]},"vetor_y", {[]},  % parâmetros vetoriais (t, Y) lidos
  "Tp_x", {-1}, "Tp_y", {-1},       % Instante de Pico
  "Ts_x", {-1}, "Ts_y", {-1},       % Tempo de acomodação
  "UP", {-1},                       % Percentual de ultrapassagem
  "zeta", {-1},                     % Taxa de amortecimento
  "Tr_xi", {-1}, "Tr_yi", {-1},     % Coord. tempo de subida (0.1 de Tp_y)
  "Tr_xf", {-1}, "Tr_yf", {-1},     % Coord. tempo de subida (0.9 de Tp_y)
  "Tr", {-1},                       % Tempo de subida
  "Wn", {-1},                       % Frequência natural
  "calculado", {sistema},           % Sub-struct dos parâmetros re-calculados.
  "funcao_de_transferencia", {""}   % Função de transferência (textual)
);

% Parâmetros do segundo sistema
sistemaB = struct (
  "nome_do_arquivo", {"empty"},     % Nome do arquivo carregado na memória
  "vetor_t", {[]},"vetor_y", {[]},  % parâmetros vetoriais (t, Y) lidos
  "Tp_x", {-1}, "Tp_y", {-1},       % Instante de Pico
  "Ts_x", {-1}, "Ts_y", {-1},       % Tempo de acomodação
  "UP", {-1},                       % Percentual de ultrapassagem
  "zeta", {-1},                     % Taxa de amortecimento
  "Tr_xi", {-1}, "Tr_yi", {-1},     % Coord. tempo de subida (0.1 de Tp_y)
  "Tr_xf", {-1}, "Tr_yf", {-1},     % Coord. tempo de subida (0.9 de Tp_y)
  "Tr", {-1},                       % Tempo de subida
  "Wn", {-1},                       % Frequência natural
  "calculado", {sistema},           % Sub-struct dos parâmetros re-calculados.
  "funcao_de_transferencia", {""}   % Função de transferência (textual)
);
                  
% FUNÇÕES E PROCEDIMENTOS ÚTEIS ------------------------------------------------

% Exibe o conteúdo de texto da variável "mensagem" na Janela do Console.
function exibirMensagem(mensagem, pularLinha = true)
  if pularLinha == true
    printf ("\a%s\n", mensagem);
  elseif
    printf ("\a%s", mensagem);
  endif
endfunction

% Verifica arquivos presentes no diretório atual
function verificarArquivos = listarArquivos()
  tmpIndex = 0;
  tmpNomes = {};
  for u = 1:15 
    nome = sprintf('Sist_2Ord%d.mat', u);
    if exist(nome, "file")
      exibirMensagem(sprintf("%s encontrado, ", nome), false);
      tmpIndex = tmpIndex + 1;
      tmpNomes{tmpIndex} = nome;
    elseif
      exibirMensagem(sprintf("%s não encontrado, ", nome), false);
    endif  
  endfor  
  verificarArquivos = tmpNomes; 
  exibirMensagem(""); % quebra de linha no console  
endfunction

% Retorna o último valor do vetor de dados fornecido
function calcValorFinalDoVetor = vetFinal(vetor)
  calcValorFinalDoVetor = vetor(end);
endfunction

% Cálculo do Instante de Pico (Tp)
function calcInstanteDePico = Tp(t, y)
  y_pmax = max(y);
  x_paux = find(y_pmax == y);
  x_pmax = t(x_paux);
  calcInstanteDePico = {x_pmax y_pmax};
endfunction

% Cálculo do tempo de acomodação (Ts) com base na tolerância (padrão 2%).
function calcTempoDeAcomodacao = Ts(t, y, tolerancia)
  y_final = y(end);
  for k=length(y):-1:0
    if((y(k)<(1-(tolerancia/100))*y_final)||(y(k)>(1+(tolerancia/100))*y_final))
      x_acomodacao = t(k-1);
      y_acomodacao = y(k-1);
      break;
    end
  end
  calcTempoDeAcomodacao = {x_acomodacao y_acomodacao};
endfunction

% Cálculo do percentual de ultrapassagem (UP)
function calcPercentualDeUltrapassagem = UP(y)
  y_final = y(end);
  y_max = max(y);
  calcPercentualDeUltrapassagem = ((y_max - y_final) / (y_final)) * 100;  
endfunction

% Cálculo da taxa de amortecimento (Zeta)
function calcZeta = Zeta(up)
  calcZeta = sqrt((log(up/100))^2/((log(up/100))^2 + (pi)^2));
endfunction

% Cálculo do tempo de subida (0.1 a 0.9), retorna coordenadas.
function calcTempoDeSubida = Tr(t, y)
  y_final = y(end);
  y_rini = 0.1 * y_final;
  y_rinia = abs(y - y_rini);
  x_rini = t(find(y_rinia == min(y_rinia)));
  y_rend = 0.9 * y_final;
  y_renda = abs(y - y_rend);
  y_rendm = find(diff(y_renda) > 0, 1, 'first');
  x_rend = t(y_rendm);
  tr = x_rend - x_rini;
  calcTempoDeSubida = {x_rini y_rini x_rend y_rend tr};
endfunction

% Cálculo da frequência natural (wn)
function calcFrequenciaNatural = Wn(tpx, zeta)
  calcFrequenciaNatural = pi/(tpx * sqrt(1-(zeta^2)));
endfunction

% Cálculo do tempo de subida normalizado.
function calcTempoDeSubidaNormalizado = Tr_norm(zeta)
  % Procedimento 1 - gera resultados satisfatóps
  aux_zeta = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
  aux_tr_norm_vet = [1.104 1.203 1.321 1.463 1.638 1.854 2.126 2.467 2.883];
  aux_zeta_a = abs(aux_zeta - zeta);
  aux_zeta_m = find(diff(aux_zeta_a)>0,1,'first'); 
  zeta
  aux_tr_norm_vet(aux_zeta_m)  
  calcTempoDeSubidaNormalizado = aux_tr_norm_vet(aux_zeta_m);
  
  % Procedimento 2 - por interpolação tem-se melhores resultados.
  % 
  % aux_zeta = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
  % aux_tr_norm_vet = [1.104 1.203 1.321 1.463 1.638 1.854 2.126 2.467 2.883];
  % p = polyfit(aux_zeta, aux_tr_norm_vet, 3);
  % ou p = (1.768*zeta^3 - 0.417 * zeta^2 + 1.039 * zeta + 1) / wn.)
  % calcTempoDeSubidaNormalizado = polyval(p, zeta);
endfunction

% CÓDIGO PRINCIPAL -------------------------------------------------------------

% Mensagem de Boas-Vindas no Console!
exibirMensagem("Procurando arquivos do tipo Sist_2Ord*.mat...");

% Verifica a presença dos arquivos Sist_2Ord*.mat no diretorio atual.
fileNames = listarArquivos();

% Continua, caso tenha encontrado arquivos, do contrário avisa o usuário
if length(fileNames) <= 0
  msgbox ("\tInfelizmente não há arquivos do tipo '.mat' neste diretório.\
  Experimente adicionar os arquivos 'Sist_2Ord*.mat' para continuar.",
  "Finalizando Programa");
elseif
  % Escolhe dois arquivos aleatoriamente e carrega suas informações preliminares
  sistemaA.nome_do_arquivo = fileNames{floor(1 + length(fileNames) * rand(1))};
  load(sistemaA.nome_do_arquivo);
  sistemaA.vetor_t = t;
  sistemaA.vetor_y = y;  
  
  sistemaB.nome_do_arquivo = fileNames{floor(1 + length(fileNames) * rand(1))};  
  load(sistemaB.nome_do_arquivo);
  sistemaB.vetor_t = t;
  sistemaB.vetor_y = y;
  
  % Limpa variáveis não utilizadas da memória
  clear -variables "t" "y";
  
  % Instante de pico (Tp)
  _tp = Tp(sistemaA.vetor_t, sistemaA.vetor_y);
  sistemaA.Tp_x = _tp{1};
  sistemaA.Tp_y = _tp{2};
  
  _tp = Tp(sistemaB.vetor_t, sistemaB.vetor_y);
  sistemaB.Tp_x = _tp{1};
  sistemaB.Tp_y = _tp{2};
  clear -variables "_tp";
  
  % Tempo de acomodação (Ts)
  _ts = Ts(sistemaA.vetor_t, sistemaA.vetor_y, TOLERANCIA_TEMPO_DE_ACOMODACAO);
  sistemaA.Ts_x = _ts{1};
  sistemaA.Ts_y = _ts{2};
  
  _ts = Ts(sistemaB.vetor_t, sistemaB.vetor_y, TOLERANCIA_TEMPO_DE_ACOMODACAO);
  sistemaB.Ts_x = _ts{1};
  sistemaB.Ts_y = _ts{2};
  clear -variables "_ts";
  
  % Percentual de Ultrapassagem
  sistemaA.UP = UP(sistemaA.vetor_y); 
  sistemaB.UP = UP(sistemaB.vetor_y);
  
  % Taxa de amortecimento (0 <= zeta < 1)
  sistemaA.zeta = Zeta(sistemaA.UP);
  sistemaB.zeta = Zeta(sistemaB.UP);
  
  % Tempo de subida (para subir de 10% à 90% da Tp_y)
  _tr = Tr(sistemaA.vetor_t, sistemaA.vetor_y);
  sistemaA.Tr_xi = _tr{1};
  sistemaA.Tr_yi = _tr{2};
  sistemaA.Tr_xf = _tr{3};
  sistemaA.Tr_yf = _tr{4};
  sistemaA.Tr = _tr{5};
  
  _tr = Tr(sistemaB.vetor_t, sistemaB.vetor_y);
  sistemaB.Tr_xi = _tr{1};
  sistemaB.Tr_yi = _tr{2};
  sistemaB.Tr_xf = _tr{3};
  sistemaB.Tr_yf = _tr{4};
  sistemaB.Tr = _tr{5};
  clear -variables "_tr";
  
  % Frequência natural (wn)
  sistemaA.Wn = Wn(sistemaA.Tp_x, sistemaA.zeta);
  sistemaB.Wn = Wn(sistemaB.Tp_x, sistemaB.zeta);
  
  % Encontra ambas as funções de transferências
  s = tf('s');
  k_dc = sistemaA.vetor_y(end);
  G1 = (k_dc * (sistemaA.Wn)^2)/(s^2 + (2*sistemaA.zeta*sistemaA.Wn*s) + 
    (sistemaA.Wn)^2);  
  sistemaA.funcao_de_transferencia = sprintf(
    'G(s)_1 = %.2f/(s^2 + %.2fs + %.2f).',
    k_dc * (sistemaA.Wn^2), 2 * sistemaA.zeta * sistemaA.Wn, sistemaA.Wn^2);
  k_dc = sistemaB.vetor_y(end);
  G2 = (k_dc * (sistemaB.Wn)^2)/(s^2 + (2*sistemaB.zeta*sistemaB.Wn*s) + 
    (sistemaB.Wn)^2);
  sistemaB.funcao_de_transferencia = sprintf(
    'G(s)_2 = %.2f/(s^2 + %.2fs + %.2f).',
    k_dc * (sistemaB.Wn^2), 2 * sistemaB.zeta * sistemaB.Wn, sistemaB.Wn^2);    
  clear -variables "k_dc";
  
  % De posse das funções de transferência, calcula-se novamente os parâmetros
  % para validação das expressões encontradas, G(s)_1 e G(S)_2, respectivamente.    
  sistemaA.calculado.Wn = sistemaA.Wn;
  sistemaA.calculado.zeta = sistemaA.zeta;
  sistemaA.calculado.Tp = pi/(sistemaA.calculado.Wn*sqrt(
    1-(sistemaA.calculado.zeta)^2));
  sistemaA.calculado.Ts = (4)/(sistemaA.calculado.zeta*sistemaA.calculado.Wn);
  sistemaA.calculado.Tr = Tr_norm(sistemaA.calculado.zeta)/(
    sistemaA.calculado.Wn);  
  sistemaA.calculado.UP = exp(-((sistemaA.calculado.zeta*pi)/(
    sqrt(1-(sistemaA.calculado.zeta)^2))))*100;
  
  sistemaB.calculado.Wn = sistemaB.Wn;
  sistemaB.calculado.zeta = sistemaB.zeta;
  sistemaB.calculado.Tp = pi/(sistemaB.calculado.Wn*sqrt(
    1-(sistemaB.calculado.zeta)^2));
  sistemaB.calculado.Ts = 4/(sistemaB.calculado.zeta*sistemaB.calculado.Wn);
  sistemaB.calculado.Tr = Tr_norm(sistemaB.calculado.zeta)/(
    sistemaB.calculado.Wn);
  sistemaB.calculado.UP = exp(-((sistemaB.calculado.zeta*pi)/(
    sqrt(1-(sistemaB.calculado.zeta)^2))))*100;
  
  % Imprime resultados no console
  exibirMensagem("");
  exibirMensagem("-----------------------------------------------------------");
  exibirMensagem(sprintf("Arquivo 1: %s", sistemaA.nome_do_arquivo));
  exibirMensagem("");exibirMensagem("Função de Transferência:");
  exibirMensagem(sprintf("%s", sistemaA.funcao_de_transferencia));
  exibirMensagem("");exibirMensagem("Parâmetros:");
  exibirMensagem(sprintf('Tp = %.2f s', sistemaA.Tp_x));
  exibirMensagem(sprintf("Ts = %.2f s", sistemaA.Ts_x));
  exibirMensagem(sprintf("Tr = %.2fs - %.2fs = %.2f s", 
    sistemaA.Tr_xf, sistemaA.Tr_xi, sistemaA.Tr));
  exibirMensagem(sprintf("Wn =  %.2f rad/s", sistemaA.Wn));
  exibirMensagem(sprintf("Zeta =  %.2f ", sistemaA.zeta));  
  exibirMensagem(sprintf("UP%% = %.2f%%", sistemaA.UP));  
  exibirMensagem("");exibirMensagem("Parâmetros Re-Calculados:");
  exibirMensagem(sprintf('Tp = %.2f s', sistemaA.calculado.Tp));
  exibirMensagem(sprintf("Ts = %.2f s", sistemaA.calculado.Ts));
  exibirMensagem(sprintf("Tr = %.2f s", sistemaA.calculado.Tr));
  exibirMensagem(sprintf("Wn =  %.2f rad/s", sistemaA.calculado.Wn));
  exibirMensagem(sprintf("Zeta =  %.2f ", sistemaA.calculado.zeta));  
  exibirMensagem(sprintf("UP%% = %.2f%%", sistemaA.calculado.UP)); 
  exibirMensagem("-----------------------------------------------------------");
  exibirMensagem("");
  exibirMensagem("-----------------------------------------------------------");
  exibirMensagem(sprintf("Arquivo 2: %s", sistemaB.nome_do_arquivo));
  exibirMensagem("");exibirMensagem("Função de Transferência:");
  exibirMensagem(sprintf("%s", sistemaB.funcao_de_transferencia));
  exibirMensagem("");exibirMensagem("Parâmetros:");
  exibirMensagem(sprintf('Tp = %.2f s', sistemaB.Tp_x));
  exibirMensagem(sprintf("Ts = %.2f s", sistemaB.Ts_x));
  exibirMensagem(sprintf("Tr = %.2fs - %.2fs = %.2f s", 
    sistemaB.Tr_xf, sistemaB.Tr_xi, sistemaB.Tr));
  exibirMensagem(sprintf("Wn =  %.2f rad/s", sistemaB.Wn));
  exibirMensagem(sprintf("Zeta =  %.2f ", sistemaB.zeta));  
  exibirMensagem(sprintf("UP%% = %.2f%%", sistemaB.UP));  
  exibirMensagem("");exibirMensagem("Parâmetros Re-Calculados:");
  exibirMensagem(sprintf('Tp = %.2f s', sistemaB.calculado.Tp));
  exibirMensagem(sprintf("Ts = %.2f s", sistemaB.calculado.Ts));
  exibirMensagem(sprintf("Tr = %.2f s", sistemaB.calculado.Tr));
  exibirMensagem(sprintf("Wn =  %.2f rad/s", sistemaB.calculado.Wn));
  exibirMensagem(sprintf("Zeta =  %.2f ", sistemaB.calculado.zeta));  
  exibirMensagem(sprintf("UP%% = %.2f%%", sistemaB.calculado.UP)); 
  exibirMensagem("-----------------------------------------------------------");
  
  % Exibe os gráficos resultantes  
  figure(1,"position",get(0,"screensize")); % Abre janela maximizada.
  subplot(2, 3, 1); plot(sistemaA.vetor_t, sistemaA.vetor_y, "-r"), hold on;
  axis([0 20 0 sistemaA.Tp_y*1.1]);
  title(strcat('\color{red}',sistemaA.nome_do_arquivo));
  plot([0 sistemaA.Tp_x],[sistemaA.Tp_y sistemaA.Tp_y],'--k','Linewidth',1.4);
  plot([sistemaA.Tp_x sistemaA.Tp_x],[0 sistemaA.Tp_y],'--k','Linewidth',1.4);
  plot([0 sistemaA.Ts_x],[sistemaA.Ts_y sistemaA.Ts_y],'--k','Linewidth',1.4);
  plot([sistemaA.Ts_x sistemaA.Ts_x],[0 sistemaA.Ts_y],'--k','Linewidth',1.4);
  plot([0 sistemaA.Tr_xi],[sistemaA.Tr_yi sistemaA.Tr_yi],
    '--k','Linewidth',1.4);
  plot([sistemaA.Tr_xi sistemaA.Tr_xi],[0 sistemaA.Tr_yi],
    '--k','Linewidth',1.4);
  plot([0 sistemaA.Tr_xf],[sistemaA.Tr_yf sistemaA.Tr_yf],
    '--k','Linewidth',1.4);
  plot([sistemaA.Tr_xf sistemaA.Tr_xf],[0 sistemaA.Tr_yf],
    '--k','Linewidth',1.4);
  text (sistemaA.Tp_x + 0.5, sistemaA.Tp_y, 
    sprintf('Tp = %.2f s', sistemaA.Tp_x));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(5+((2*5)/5)), 
    sprintf("Ts = %.2f s", sistemaA.Ts_x));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(4+((2*5)/5)), 
    sprintf("Tr = %.2fs - %.2fs = %.2f s", 
      sistemaA.Tr_xf, sistemaA.Tr_xi, sistemaA.Tr));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(3+((2*5)/5)), 
    strcat('\omega_n = ',sprintf("  %.2f rad/s", sistemaA.Wn)));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(2+((2*5)/5)), 
    strcat('\zeta = ',sprintf(" %.2f ", sistemaA.zeta)));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(1+((2*5)/5)), 
    sprintf("UP%% = %.2f %% ", sistemaA.UP));  
  legend("off"); xlabel("t"); ylabel("y");grid off;
  subplot(2, 3, 2); step(G1, "-b", [0:0.01:20]), hold on;
  axis([0 20 0 sistemaA.Tp_y*1.1]);
  text (sistemaA.Tp_x + 0.5, sistemaA.Tp_y, 
    sprintf("Tp = %.2f s", sistemaA.calculado.Tp));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(5+((2*5)/5)), 
    sprintf("Ts = %.2f s", sistemaA.calculado.Ts));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(4+((2*5)/5)), 
    sprintf("Tr = %.2f s", sistemaA.calculado.Tr));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(3+((2*5)/5)), 
    strcat('\omega_n = ', sprintf(" %.2f rad/s", sistemaA.calculado.Wn)));
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(2+((2*5)/5)), 
    strcat('\zeta = ',sprintf(" %.2f ", sistemaA.calculado.zeta)));  
  text (sistemaA.Ts_x + 0.5, (sistemaA.Ts_y/(2*5))*(1+((2*5)/5)), 
    sprintf("UP%% = %.2f %% ", sistemaA.calculado.UP));  
  title(strcat('\color{blue}',sistemaA.funcao_de_transferencia));
  legend("off");xlabel("t");ylabel("y");grid off;
  subplot(2, 3, 3); plot(sistemaA.vetor_t, sistemaA.vetor_y, "-g"), hold on;    
  step(G1, "-xb", [0:0.5:20]), hold on;
  axis([0 20 0 sistemaA.Tp_y*1.1]);
  plot(sistemaA.vetor_t, sistemaA.vetor_y, "-r"), hold on;  
  title(strcat('\color{red}',sistemaA.nome_do_arquivo,'\color{black} versus ',
    '\color{blue} G(s)_1'));
  legend("off");xlabel("t");ylabel("y");grid off;
  
  subplot(2, 3, 4); plot(sistemaB.vetor_t, sistemaB.vetor_y, "-r"), hold on;  
  axis([0 20 0 sistemaB.Tp_y*1.1]);
  title(strcat('\color{red}', sistemaB.nome_do_arquivo));xlabel("Tempo (s)");
  plot([0 sistemaB.Tp_x],[sistemaB.Tp_y sistemaB.Tp_y],'--k','Linewidth',1.4);
  plot([sistemaB.Tp_x sistemaB.Tp_x],[0 sistemaB.Tp_y],'--k','Linewidth',1.4);
  plot([0 sistemaB.Ts_x],[sistemaB.Ts_y sistemaB.Ts_y],'--k','Linewidth',1.4);
  plot([sistemaB.Ts_x sistemaB.Ts_x],[0 sistemaB.Ts_y],'--k','Linewidth',1.4);
  plot([0 sistemaB.Tr_xi],[sistemaB.Tr_yi sistemaB.Tr_yi],
    '--k','Linewidth',1.4);
  plot([sistemaB.Tr_xi sistemaB.Tr_xi],[0 sistemaB.Tr_yi],
    '--k','Linewidth',1.4);
  plot([0 sistemaB.Tr_xf],[sistemaB.Tr_yf sistemaB.Tr_yf],
    '--k','Linewidth',1.4);
  plot([sistemaB.Tr_xf sistemaB.Tr_xf],[0 sistemaB.Tr_yf],
    '--k','Linewidth',1.4);
  text (sistemaB.Tp_x + 0.5, sistemaB.Tp_y, 
    sprintf('Tp = %.2f s', sistemaB.Tp_x));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(5+((2*5)/5)), 
    sprintf("Ts = %.2f s", sistemaB.Ts_x));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(4+((2*5)/5)), 
    sprintf("Tr = %.2fs - %.2fs = %.2f s", 
      sistemaB.Tr_xf, sistemaB.Tr_xi, sistemaB.Tr));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(3+((2*5)/5)), 
    strcat('\omega_n = ',sprintf(' %.2f rad/s', sistemaB.Wn)));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(2+((2*5)/5)), 
    strcat('\zeta = ', sprintf(" %.2f ", sistemaB.zeta)));  
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(1+((2*5)/5)), 
    sprintf("UP%% = %.2f %% ", sistemaB.UP));  
  legend("off"); xlabel("t");ylabel("y");grid off;
  subplot(2, 3, 5); step(G2, "-b", [0:0.01:20]), hold on;
  axis([0 20 0 sistemaB.Tp_y*1.1]);
  text (sistemaB.Tp_x + 0.5, sistemaB.Tp_y, sprintf("Tp = %.2f s", 
    sistemaB.calculado.Tp));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(5+((2*5)/5)), 
    sprintf("Ts = %.2f s", sistemaB.calculado.Ts));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(4+((2*5)/5)), 
    sprintf("Tr = %.2f s", sistemaB.calculado.Tr));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(3+((2*5)/5)), 
    strcat('\omega_n = ', sprintf(" %.2f rad/s", sistemaB.calculado.Wn)));
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(2+((2*5)/5)), 
    strcat('\zeta = ',sprintf(" %.2f ", sistemaB.calculado.zeta)));  
  text (sistemaB.Ts_x + 0.5, (sistemaB.Ts_y/(2*5))*(1+((2*5)/5)), 
    sprintf("UP%% = %.2f %% ", sistemaB.calculado.UP));  
  title(strcat('\color{blue}', sistemaB.funcao_de_transferencia));
  legend("off");xlabel("t");ylabel("y");grid off;
  subplot(2, 3, 6); step(G2, "-xb", [0:0.5:20]), hold on;
  axis([0 20 0 sistemaB.Tp_y*1.1]);
  plot(sistemaB.vetor_t, sistemaB.vetor_y, "-r"), hold on;  
  title(strcat('\color{red}',sistemaB.nome_do_arquivo,'\color{black} versus ',
    '\color{blue} G(s)_2'));
  legend("off");xlabel("t");ylabel("y");grid off;
  %text (1, 1, 'exemplo');
  %annotation('Textarrow',[0.7 0.7],[0.7 0.52],'FontSize',13,'Linewidth',2);
  hold off;
endif
