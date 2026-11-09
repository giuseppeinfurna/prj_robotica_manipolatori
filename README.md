# Robot manipolatore planare 4R con gravità

Questo progetto MATLAB modella un braccio robotico planare a quattro giunti rotoidali che si muove nel piano XY, con la gravità diretta lungo -y. Partendo dalla cinematica diretta e inversa, passa per la dinamica ricavata con Eulero-Lagrange e arriva a un controllore a coppia calcolata che fa seguire all'end-effector una traiettoria circolare. Le unità sono quelle del Sistema Internazionale, gli angoli sono in radianti e le rotazioni positive sono antiorarie attorno all'asse z uscente dal piano.

## Come eseguirlo

Basta mettere tutti i file `.m` nella stessa cartella, impostarla come cartella corrente in MATLAB e lanciare `main_4R`. Lo script verifica prima la cinematica e la dinamica, poi simula il robot con `ode45` e infine mostra i grafici e l'animazione. Per la sola derivazione simbolica serve il Symbolic Math Toolbox: si lancia `derive_dynamics_symbolic`, che ricava `M`, `C` e `G` e li confronta con la versione numerica.

Da terminale si può usare `matlab -batch "main_4R"`. In questo modo però le figure non restano aperte; per vederle conviene usare l'interfaccia grafica oppure `matlab -r "main_4R"`.

## Come è organizzato il codice

La cinematica diretta è in `fk_4R`, che moltiplica le quattro matrici di Denavit-Hartenberg standard e restituisce posizione e orientamento dell'end-effector insieme a tutte le terne intermedie. Il robot ha quattro giunti ma cinque terne: la terna 0 è fissa alla base, mentre la terna `i` è quella associata all'estremo del link `i`, quindi la quarta coincide con l'end-effector.

La cinematica inversa è in `ik_4R`. Con tre coordinate di task (x, y e l'orientamento phi) e quattro giunti il robot è ridondante di un grado, perciò occorre fissare una scelta: qui si assegna il valore di `theta4`, poi si risale al polso, si risolve il problema 2R con la legge dei coseni e si ottengono `theta1`, `theta2` e `theta3`. Le due soluzioni di gomito (alto e basso) si scelgono con un'opzione, e se il punto non è raggiungibile la funzione solleva l'errore `ik_4R:Unreachable`, con un messaggio che indica la posa richiesta.

La dinamica è in `dynamics_4R`, che costruisce la matrice di inerzia a partire dai Jacobiani dei baricentri, il vettore di gravità dal gradiente dell'energia potenziale e i termini di Coriolis e centrifughi con i simboli di Christoffel. Gli stessi termini sono ricavati in forma simbolica da `derive_dynamics_symbolic`, che serve da controprova della versione numerica.

Il disegno del robot con le terne DH è affidato a `plot_robot_4R`, usato sia per le figure statiche sia per l'animazione in `main_4R`.

## Verifiche automatiche

Prima di simulare, `main_4R` controlla con `assert` che il modello sia coerente, e si ferma con un messaggio chiaro se qualcosa non torna. Verifica che `phi` sia la somma dei giunti e che la parte rotazionale di `T` sia ortogonale, che `FK(IK(p))` riproduca il punto di partenza per entrambi i gomiti, che un punto fuori dallo spazio di lavoro generi l'errore previsto, che `M` sia simmetrica e definita positiva e che l'energia si conservi nel moto libero senza coppie. Al termine della simulazione controlla anche che l'errore di inseguimento, dopo il transitorio iniziale, resti sotto un centesimo di radiante. La derivazione simbolica fa lo stesso confronto tra `M`, `C` e `G` simbolici e numerici, e verifica anche che `Mdot - 2C` sia antisimmetrica.

La figura con `q(t)`, l'errore, le coppie e la traiettoria dell'end-effector viene salvata automaticamente in `results/risultati_controllo.png` a ogni esecuzione.

## Risultati della validazione

Il codice è stato eseguito su MATLAB Online (R2026b) e tutti i controlli preliminari sono stati superati. Per la configurazione di test, con giunti a 30, -40, 25 e 15 gradi, l'end-effector si trova in x = 1.2899 m e y = 0.3582 m con orientamento phi = 0.5236 rad. La cinematica inversa ritrova la stessa posa con entrambi i gomiti: il gomito alto restituisce esattamente la configurazione di partenza, mentre quello basso ne trova una diversa (-5.37, 40.00, -19.63 e 15.00 gradi) che raggiunge lo stesso punto. L'errore di ricostruzione FK(IK(p)) è dell'ordine di 1e-16, cioè al limite della precisione di macchina.

Nel moto libero senza coppie l'energia totale, che parte da 6.82 J, varia al massimo di circa 1.2e-08 J, valore compatibile con la tolleranza di `ode45` e quindi coerente con un modello dinamico corretto.

Con il controllore a coppia calcolata, partendo da un errore di giunto iniziale di circa 0.25 rad, il robot converge sulla circonferenza entro due secondi: dopo quel momento l'errore resta sotto 1.2e-08 rad, e a fine simulazione scende a 1.6e-11 rad sui giunti e a 5.5e-12 m sulla posizione dell'end-effector. Le coppie massime richieste sono di circa 13.4, 11.9, 3.5 e 0.5 N m per i quattro giunti, dalla base verso l'estremo.

## Controllore e traiettoria

Il riferimento è una circonferenza nello spazio operativo, trasformata in traiettoria di giunto tramite la cinematica inversa e interpolata con una spline per ottenere velocità e accelerazioni. La legge di controllo è `tau = M(q)(ddq_d + Kd e' + Kp e) + C(q,dq) + G(q)`, che rende l'errore un sistema lineare asintoticamente stabile; i guadagni sono scelti per avere smorzamento critico. Il modello usato dal controllore coincide con quello del robot simulato, quindi non sono considerati attriti, saturazioni delle coppie o errori di modello.

## Una versione alternativa

Il file `sim_4R.m` contiene la stessa idea in un unico script, nello stile delle esercitazioni sincrone del corso: una variabile `CONTROL_MODE` permette di scegliere tra controllo nello spazio dei giunti e nello spazio operativo, e la dinamica è integrata con il metodo di Eulero invece che con `ode45`. Serve soprattutto come confronto didattico; per la consegna fa fede `main_4R`.
