clear;
% Name: Hamza Malik, Student number: 501112545, Section: 12
%%
% ELE532_LAB_1: Working with Matlab Functions, Visualization of Signals, and Signals Properties
%
%% Part A: Anonymous functions and plotting continuous functions
%Problem A.1: 

    %Figure 1.46
t1 = (-2:2);
f1= @(t1) exp(-t1).*cos(2*pi*t1);
figure;
plot(t1,f1(t1));
xlabel('time'); ylabel('f1(t1)'); title('Figure 1.46'); grid;
hold off;

    %Figure 1.47
t2 = (-2:0.01:2);
f2= @(t2) exp(-t2).*cos(2*pi*t2);
figure;
plot(t2,f2(t2));
xlabel('time'); ylabel('f2(t2)'); title('Figure 1.47'); grid;
hold off;

%Problem A.2: 
t = (-2:1:2);
y = exp(-t);
figure;
plot(t,y);
xlabel('time'); ylabel('y(t)'); title('Function A.2'); grid;
hold off;

%Problem A.3:

    %The function plotted in problem A.2 is the exact same as the visual 
    %representation of the function of figure 1.46 in problem A.1

%% Part B: Time shifting and time scaling
%Problem B.1
t = (-1:0.01:2);
p = @(t) 1.0.*((t>=0)&(t<1));
figure;
plot(t,p(t));
xlabel('time'); ylabel('p(t)'); title('Figure 1.50'); grid;
hold off;

%Problem B.2
r = @(t) p(t) .* t;
figure;
plot(t,r(t));
xlabel('time'); ylabel('r(t)'); title('r(t) PLOT'); grid;
hold off;

n = @(t) r(t) + r(-t+2);
figure;
plot(t,n(t));
xlabel('time'); ylabel('n(t)'); title('n(t) PLOT'); grid;
hold off;

%Problem B.3
n1 = @(t) t .* 0.5;
n2 = @(t) n1(t + 0.5);
figure;
plot(t, n1(t), t, n2(t));
xlabel('time'); ylabel('n(t)'); title('n1(t) & n2(t)'); grid;
hold on;

%Problem B.4 
n3 = @(t) t + 0.25;
n4 = @(t) n3(t .* 0.5);
figure;
plot(t, n3(t), t, n4(t));
xlabel('time'); ylabel('n(t)'); title('n3(t) & n4(t)'); grid;
hold on;

%Problem B.5
figure;
plot(t, n2(t), t, n4(t));
xlabel('time'); ylabel('n(t)'); title('n2(t) & n4(t)'); grid;
hold on;
    %Graphs n2(t) and n4(t) are the exact same graph, sharing the same
    %slopes and intersects.


%% Part C: Time shifting and time scaling

%Problem C.1
t = (-2:0.01:2);
f = @(t) exp(-2.*t).*cos(4.*pi.*t);
u = @(t) 1.0.*(t>=0);
g = @(t) f(t).*u(t);
figure;
plot(t, g(t), t, f(t));
xlabel('time'); ylabel('g(t)'); title('g(t) & f(t)'); grid;
hold on;

%Problem C.2
T = (-2:0.01:4);
p = @(T) exp(-2).*g(T+1);
figure;
plot(T, p(T));
xlabel('T'); ylabel('p(T)'); title('Graph of p(T)'); grid;
hold off;

% Problem C.3
u = @(t) 1.0.*(t>=0);
t = (1:0.01:4); 
a = [1, 3, 5, 7];
s_alpha = zeros(length(a), length(t));

for alpha = a
    s = @(t) exp(-2).*exp(-alpha.*t).*cos(4*pi*t).*u(t);
    plot(t,s(t));
    hold on;
end
xlabel("t"); ylabel("a(t)"); title("Graph of s_a(t)"); grid;
legend('alpha = 1', 'alpha = 3', 'alpha = 5', 'alpha = 7');
hold off;

%Problem C.4
t = (1:0.01:4); 
alpha_matrix = a; 
t_matrix = t;
% Matrix would be 4 x 399 if transposed properly

%% Part D: Array Indexing

A = [0.5377, -1.3077, -1.3499, -0.2050; 
     1.8339, -0.4336, 3.0349, -0.1241;
     -2.2588 0.3426 0.7254 1.4897;
     0.8622 3.5784, -0.0631 1.4090;
     0.3188 2.7694 0.7147 1.4172;];

%Problem D.1
A(:)
A([2 4 7])
[A>=0.2]
A([A >= 0.2])
A([A>=0.2]) = 0;

% ✅ Fix for Problem D.2: define B if it doesn't exist
if ~exist('B', 'var')
    B = randn(1024, 100);  % Replace with actual data if available
end

Num_rows = size(B,1);
Num_cols = size(B,2);

for i = 1:Num_rows
    for j = 1:Num_cols
        if abs(B(i,j)) < 0.01
            B(i,j) = 0;
        end
    end
end

B([abs(B)>= 0.01]) = 0;

% part Ci
tic
for i = 1:Num_rows
    for j = 1:Num_cols
        if abs(B(i,j)) < 0.01
            B(i,j) = 0;
        end
    end
end
fprintf('\nPART A execution time: ')
toc

% part Cii
tic
B([abs(B)>= 0.01]) = 0;
fprintf('PART B execution time: ')
toc

% Problem D.3
if ~exist('x_audio', 'var')
    x_audio = randn(20000,1);  % Replace with actual audio data if available
end

Num_rows = size(x_audio,1);
Num_cols = size(x_audio,2);
threshold = 0;

for i = 1:Num_rows
    for j = 1:Num_cols
        if abs(x_audio(i,j)) == 0
            threshold = threshold + 1;
        end
    end
end

f = @(x_audio) sum(~x_audio(:));
f(x_audio)
fprintf("For the audio data set, the threshold is: " + threshold);

sound(x_audio,8000);


