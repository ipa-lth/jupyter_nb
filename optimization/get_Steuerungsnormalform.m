% https://github.com/ipa-lth/jupyter_nb/blob/master/optimization/optim_tools.py
function [A0, b0, c0, d, T, Q] = get_Steuerungsnormalform(A, b, c, d)  
  #disp ("Running get_Steuerungsnormalform...");
  #https://www.eit.hs-karlsruhe.de/mesysto/teil-a-zeitkontinuierliche-signale-und-systeme/darstellung-von-systemen-im-zustandsraum/transformation-auf-eine-bestimmte-darstellungsform/transformation-einer-zustandsgleichung-in-regelungsnormalform.html

  # Image(url = "https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_7_HQ.png")


  # Berechnung der inversen Steuerbarkeitsmatrix
  n = size (A, 1);
  #disp ("n = ");
  #disp(n);
  Q = b;
  for i = 1: (n-1)
    Q = [Q, (A^i)*b];
  endfor
  #disp ("Q = ");
  #disp(Q);
  Q_inv = (Q)^-1;
  #disp ("Q_inv = ");
  #disp(Q_inv);

  #Zeilenvektor t_1.T entspricht der letzten Zeile der inversen Steuerbarkeitsmatrix
  #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Formel_10_3_51_HQ.png")

  t1 = Q_inv(size (Q_inv, 1),:);

  # Berechnung der Transformationsmatrix 
  #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_8_HQ.png")
  T = t1;
  for i = 1: (n-1)
      T = [T; t1*(A^i)];
  endfor
  #Bestimmung der Zustandsraumdarstellung in Regelungsnormalform 
  #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_9_HQ.png")
  #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_10_HQ.png")

  A0 = T*A*(T^-1);
  b0 = T*b;
  c0 = (c.' * (T^-1)).';
endfunction