%% EGB342 Assignment 2B - Data Generator
%  Use this file to generate the data that you will be using in
%  Assignment 2B.
%  You will only need to run this code ONCE.
%  Three files will be generated, which includes
%  all of the data that you will be using in this assignment.
%  - Fill in your group number on line 15.
%  - Fill in student numbers on lines 24-26, as per instructions.
%  - Then run the code.
%
%% Preparing MATLAB workspace
clear all, clc, close all;

%% Student input %%
% ----- Enter you group number here ----- %
g = 7;   % <<<-------------------------------------------- Group Number
%  
%  Student numbers have the format "n01234567".
%  Omit the leading 'n' and leading '0', and enter it as 1234567.
%  It should only be a 7 digit number.
%  If there are only two people in your group, take the sum of the student
%  numbers and discard the leading digit. Then insert the result into st3.

st1 = 9630210; % Student 1   % <<<----------------------- Student Numbers
st2 = 9698701; % Student 2 
st3 = 9713310; % Student 3

% End of student input %%



%% STOP STOP STOP

%% Do not chage the code below
%% Part 1
TextGenerator(st1,st2,st3);

fprintf('Your text file for Part 1 has been saved to:\n\r');
disp(pwd)
fprintf('In the file "sample.txt".\n\r');

%% Part 2
st = [int2str(st1),' ',int2str(st2),' ',int2str(st3)];
% Generating data
data = A2BGenData(g,st);
x = data.xdata;
qpsk_msg = data.qpsk_msg;
test_msg_str = data.test_msg_str;

save A2BPart2 x;
save A2BPart3 x qpsk_msg test_msg_str;

fprintf('Your .mat files for Part 2 and Part 3 has been saved to:\n\r');
disp(pwd)
fprintf('In the files "A2BPart2.mat and A2BPart3.mat".\n\r');
fprintf('Use the "load" command to load the saved variables into the workspace.\n');

clear all
