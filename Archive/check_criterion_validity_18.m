% This script can be used to check criterion validity for network metrics
% on the 18-term dataset

%% Instructions:
% run the following scripts for the 18-term dataset:
% - assign_subject_communities_loop (set subjectFiles to 0)
% - check_subject_communities_loop (set subjectFiles to 0)
% - check_subject_clustering_loop (set subjectFiles to 0)
% - check_subject_participation_loop (set subjectFiles to 0)
% - check_subject_density_loop
% - granularity_calculations
% - emodiversity_calculation
% - affect_calculation
% make sure the necessary .mat files are on your path

clear;
clc;

%% specify dataset
dataSet = 18; % must be set to 18
print = 0; % set to 1 to generate spreadsheet with results

%% load data and list of checks
filenames{1} = ['dataSet' num2str(dataSet) '_zICC.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_communities.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_modularity.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_clustering.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_participation.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{8} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{9} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{10} = ['dataSet' num2str(dataSet) '_comm_radius.mat'];
filenames{11} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
filenames{12} = ['dataSet' num2str(dataSet) '_affect.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end
checks = readtable('criterionChecks_18.csv');

%% load criterion validity measures
rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
ASItotal = rawData.data(:,7)-21; % subtract 21 to bring in line with original rating scale
ERStotal = rawData.data(:,11); 
GADtotal = rawData.data(:,12); 
PSStotal = rawData.data(:,21).*3.5; %46 for session 2 % multiply by 3.5 to bring in line with 14-item measure
PHQ15total = rawData.data(:,22); %47 for session 2
PHQ8total = rawData.data(:,23); %48 for session 2
neuroticism = rawData.data(:,49);


   
%% correlations between ASI and network metrics
[corr18(1,1),corr18(1,2)] = corr(ASItotal(:),numCommunities(:,2),'rows','complete');
[corr18(2,1),corr18(2,2)] = corr(ASItotal(:),numCommunities(:,3),'rows','complete');
[corr18(3,1),corr18(3,2)] = corr(ASItotal(:),numCommunities(:,4),'rows','complete');
[corr18(144,1),corr18(144,2)] = corr(ASItotal(:),modularity(:,2),'rows','complete');
[corr18(145,1),corr18(145,2)] = corr(ASItotal(:),modularity(:,3),'rows','complete');
[corr18(146,1),corr18(146,2)] = corr(ASItotal(:),modularity(:,4),'rows','complete');
[corr18(147,1),corr18(147,2)] = corr(ASItotal(:),modularity_VA(:,2),'rows','complete');
[corr18(148,1),corr18(148,2)] = corr(ASItotal(:),modularity_VA(:,3),'rows','complete');
[corr18(149,1),corr18(149,2)] = corr(ASItotal(:),modularity_VA(:,4),'rows','complete');
[corr18(4,1),corr18(4,2)] = corr(ASItotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(5,1),corr18(5,2)] = corr(ASItotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(6,1),corr18(6,2)] = corr(ASItotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(7,1),corr18(7,2)] = corr(ASItotal(:),networkDensity(:,3),'rows','complete');
[corr18(8,1),corr18(8,2)] = corr(ASItotal(:),networkDensity(:,4),'rows','complete');
[corr18(141,1),corr18(141,2)] = corr(ASItotal(:),networkDensity(:,5),'rows','complete');
[corr18(142,1),corr18(142,2)] = corr(ASItotal(:),networkDensity(:,6),'rows','complete');
[corr18(143,1),corr18(143,2)] = corr(ASItotal(:),networkDensity(:,7),'rows','complete');
[corr18(9,1),corr18(9,2)] = corr(ASItotal(:),participation(:,2),'rows','complete');
[corr18(10,1),corr18(10,2)] = corr(ASItotal(:),participation(:,3),'rows','complete');
[corr18(11,1),corr18(11,2)] = corr(ASItotal(:),participation(:,5),'rows','complete');
[corr18(12,1),corr18(12,2)] = corr(ASItotal(:),participation(:,6),'rows','complete');
[corr18(13,1),corr18(13,2)] = corr(ASItotal(:),gateway(:,2),'rows','complete'); % node strength
[corr18(14,1),corr18(14,2)] = corr(ASItotal(:),gateway(:,3),'rows','complete'); % node strength
[corr18(15,1),corr18(15,2)] = corr(ASItotal(:),gateway(:,5),'rows','complete'); % node strength
[corr18(16,1),corr18(16,2)] = corr(ASItotal(:),gateway(:,6),'rows','complete'); % node strength
[corr18(17,1),corr18(17,2)] = corr(ASItotal(:),diversity(:,2),'rows','complete');
[corr18(18,1),corr18(18,2)] = corr(ASItotal(:),diversity(:,3),'rows','complete');
[corr18(19,1),corr18(19,2)] = corr(ASItotal(:),diversity(:,5),'rows','complete');
[corr18(20,1),corr18(20,2)] = corr(ASItotal(:),diversity(:,6),'rows','complete');
[corr18(239,1),corr18(239,2)] = corr(ASItotal(:),commRadius(:,2),'rows','complete');
[corr18(240,1),corr18(240,2)] = corr(ASItotal(:),commRadius(:,3),'rows','complete');
[corr18(241,1),corr18(241,2)] = corr(ASItotal(:),commRadius(:,4),'rows','complete');

%% correlations between ERS and network metrics
[corr18(21,1),corr18(21,2)] = corr(ERStotal(:),numCommunities(:,2),'rows','complete');
[corr18(22,1),corr18(22,2)] = corr(ERStotal(:),numCommunities(:,3),'rows','complete');
[corr18(23,1),corr18(23,2)] = corr(ERStotal(:),numCommunities(:,4),'rows','complete');
[corr18(150,1),corr18(150,2)] = corr(ERStotal(:),modularity(:,2),'rows','complete');
[corr18(151,1),corr18(151,2)] = corr(ERStotal(:),modularity(:,3),'rows','complete');
[corr18(152,1),corr18(152,2)] = corr(ERStotal(:),modularity(:,4),'rows','complete');
[corr18(153,1),corr18(153,2)] = corr(ERStotal(:),modularity_VA(:,2),'rows','complete');
[corr18(154,1),corr18(154,2)] = corr(ERStotal(:),modularity_VA(:,3),'rows','complete');
[corr18(155,1),corr18(155,2)] = corr(ERStotal(:),modularity_VA(:,4),'rows','complete');
[corr18(24,1),corr18(24,2)] = corr(ERStotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(25,1),corr18(25,2)] = corr(ERStotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(26,1),corr18(26,2)] = corr(ERStotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(27,1),corr18(27,2)] = corr(ERStotal(:),networkDensity(:,3),'rows','complete');
[corr18(28,1),corr18(28,2)] = corr(ERStotal(:),networkDensity(:,4),'rows','complete');
[corr18(156,1),corr18(156,2)] = corr(ERStotal(:),networkDensity(:,5),'rows','complete');
[corr18(157,1),corr18(157,2)] = corr(ERStotal(:),networkDensity(:,6),'rows','complete');
[corr18(158,1),corr18(158,2)] = corr(ERStotal(:),networkDensity(:,7),'rows','complete');
[corr18(29,1),corr18(29,2)] = corr(ERStotal(:),participation(:,2),'rows','complete');
[corr18(30,1),corr18(30,2)] = corr(ERStotal(:),participation(:,3),'rows','complete');
[corr18(31,1),corr18(31,2)] = corr(ERStotal(:),participation(:,5),'rows','complete');
[corr18(32,1),corr18(32,2)] = corr(ERStotal(:),participation(:,6),'rows','complete');
[corr18(33,1),corr18(33,2)] = corr(ERStotal(:),gateway(:,2),'rows','complete'); % node strength
[corr18(34,1),corr18(34,2)] = corr(ERStotal(:),gateway(:,3),'rows','complete'); % node strength
[corr18(35,1),corr18(35,2)] = corr(ERStotal(:),gateway(:,5),'rows','complete'); % node strength
[corr18(36,1),corr18(36,2)] = corr(ERStotal(:),gateway(:,6),'rows','complete'); % node strength
[corr18(37,1),corr18(37,2)] = corr(ERStotal(:),diversity(:,2),'rows','complete');
[corr18(38,1),corr18(38,2)] = corr(ERStotal(:),diversity(:,3),'rows','complete');
[corr18(39,1),corr18(39,2)] = corr(ERStotal(:),diversity(:,5),'rows','complete');
[corr18(40,1),corr18(40,2)] = corr(ERStotal(:),diversity(:,6),'rows','complete');
[corr18(242,1),corr18(242,2)] = corr(ERStotal(:),commRadius(:,2),'rows','complete');
[corr18(243,1),corr18(243,2)] = corr(ERStotal(:),commRadius(:,3),'rows','complete');
[corr18(244,1),corr18(244,2)] = corr(ERStotal(:),commRadius(:,4),'rows','complete');

%% correlations between GAD and network metrics
[corr18(41,1),corr18(41,2)] = corr(GADtotal(:),numCommunities(:,2),'rows','complete');
[corr18(42,1),corr18(42,2)] = corr(GADtotal(:),numCommunities(:,3),'rows','complete');
[corr18(43,1),corr18(43,2)] = corr(GADtotal(:),numCommunities(:,4),'rows','complete');
[corr18(159,1),corr18(159,2)] = corr(GADtotal(:),modularity(:,2),'rows','complete');
[corr18(160,1),corr18(160,2)] = corr(GADtotal(:),modularity(:,3),'rows','complete');
[corr18(161,1),corr18(161,2)] = corr(GADtotal(:),modularity(:,4),'rows','complete');
[corr18(162,1),corr18(162,2)] = corr(GADtotal(:),modularity_VA(:,2),'rows','complete');
[corr18(163,1),corr18(163,2)] = corr(GADtotal(:),modularity_VA(:,3),'rows','complete');
[corr18(164,1),corr18(164,2)] = corr(GADtotal(:),modularity_VA(:,4),'rows','complete');
[corr18(44,1),corr18(44,2)] = corr(GADtotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(45,1),corr18(45,2)] = corr(GADtotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(46,1),corr18(46,2)] = corr(GADtotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(47,1),corr18(47,2)] = corr(GADtotal(:),networkDensity(:,3),'rows','complete');
[corr18(48,1),corr18(48,2)] = corr(GADtotal(:),networkDensity(:,4),'rows','complete');
[corr18(165,1),corr18(165,2)] = corr(GADtotal(:),networkDensity(:,5),'rows','complete');
[corr18(166,1),corr18(166,2)] = corr(GADtotal(:),networkDensity(:,6),'rows','complete');
[corr18(167,1),corr18(167,2)] = corr(GADtotal(:),networkDensity(:,7),'rows','complete');
[corr18(49,1),corr18(49,2)] = corr(GADtotal(:),participation(:,2),'rows','complete');
[corr18(50,1),corr18(50,2)] = corr(GADtotal(:),participation(:,3),'rows','complete');
[corr18(51,1),corr18(51,2)] = corr(GADtotal(:),participation(:,5),'rows','complete');
[corr18(52,1),corr18(52,2)] = corr(GADtotal(:),participation(:,6),'rows','complete');
[corr18(53,1),corr18(53,2)] = corr(GADtotal(:),gateway(:,2),'rows','complete'); % node strength
[corr18(54,1),corr18(54,2)] = corr(GADtotal(:),gateway(:,3),'rows','complete'); % node strength
[corr18(55,1),corr18(55,2)] = corr(GADtotal(:),gateway(:,5),'rows','complete'); % node strength
[corr18(56,1),corr18(56,2)] = corr(GADtotal(:),gateway(:,6),'rows','complete'); % node strength
[corr18(57,1),corr18(57,2)] = corr(GADtotal(:),diversity(:,2),'rows','complete');
[corr18(58,1),corr18(58,2)] = corr(GADtotal(:),diversity(:,3),'rows','complete');
[corr18(59,1),corr18(59,2)] = corr(GADtotal(:),diversity(:,5),'rows','complete');
[corr18(60,1),corr18(60,2)] = corr(GADtotal(:),diversity(:,6),'rows','complete');
[corr18(245,1),corr18(245,2)] = corr(GADtotal(:),commRadius(:,2),'rows','complete');
[corr18(246,1),corr18(246,2)] = corr(GADtotal(:),commRadius(:,3),'rows','complete');
[corr18(247,1),corr18(247,2)] = corr(GADtotal(:),commRadius(:,4),'rows','complete');

%% correlations between PSS and network metrics
[corr18(61,1),corr18(61,2)] = corr(PSStotal(:),numCommunities(:,2),'rows','complete');
[corr18(62,1),corr18(62,2)] = corr(PSStotal(:),numCommunities(:,3),'rows','complete');
[corr18(63,1),corr18(63,2)] = corr(PSStotal(:),numCommunities(:,4),'rows','complete');
[corr18(168,1),corr18(168,2)] = corr(PSStotal(:),modularity(:,2),'rows','complete');
[corr18(169,1),corr18(169,2)] = corr(PSStotal(:),modularity(:,3),'rows','complete');
[corr18(170,1),corr18(170,2)] = corr(PSStotal(:),modularity(:,4),'rows','complete');
[corr18(171,1),corr18(171,2)] = corr(PSStotal(:),modularity_VA(:,2),'rows','complete');
[corr18(172,1),corr18(172,2)] = corr(PSStotal(:),modularity_VA(:,3),'rows','complete');
[corr18(173,1),corr18(173,2)] = corr(PSStotal(:),modularity_VA(:,4),'rows','complete');
[corr18(64,1),corr18(64,2)] = corr(PSStotal(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(65,1),corr18(65,2)] = corr(PSStotal(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(66,1),corr18(66,2)] = corr(PSStotal(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(67,1),corr18(67,2)] = corr(PSStotal(:),networkDensity(:,3),'rows','complete');
[corr18(68,1),corr18(68,2)] = corr(PSStotal(:),networkDensity(:,4),'rows','complete');
[corr18(174,1),corr18(174,2)] = corr(PSStotal(:),networkDensity(:,5),'rows','complete');
[corr18(175,1),corr18(175,2)] = corr(PSStotal(:),networkDensity(:,6),'rows','complete');
[corr18(176,1),corr18(176,2)] = corr(PSStotal(:),networkDensity(:,7),'rows','complete');
[corr18(69,1),corr18(69,2)] = corr(PSStotal(:),participation(:,2),'rows','complete');
[corr18(70,1),corr18(70,2)] = corr(PSStotal(:),participation(:,3),'rows','complete');
[corr18(71,1),corr18(71,2)] = corr(PSStotal(:),participation(:,5),'rows','complete');
[corr18(72,1),corr18(72,2)] = corr(PSStotal(:),participation(:,6),'rows','complete');
[corr18(73,1),corr18(73,2)] = corr(PSStotal(:),gateway(:,2),'rows','complete'); % node strength
[corr18(74,1),corr18(74,2)] = corr(PSStotal(:),gateway(:,3),'rows','complete'); % node strength
[corr18(75,1),corr18(75,2)] = corr(PSStotal(:),gateway(:,5),'rows','complete'); % node strength
[corr18(76,1),corr18(76,2)] = corr(PSStotal(:),gateway(:,6),'rows','complete'); % node strength
[corr18(77,1),corr18(77,2)] = corr(PSStotal(:),diversity(:,2),'rows','complete');
[corr18(78,1),corr18(78,2)] = corr(PSStotal(:),diversity(:,3),'rows','complete');
[corr18(79,1),corr18(79,2)] = corr(PSStotal(:),diversity(:,5),'rows','complete');
[corr18(80,1),corr18(80,2)] = corr(PSStotal(:),diversity(:,6),'rows','complete');
[corr18(248,1),corr18(248,2)] = corr(PSStotal(:),commRadius(:,2),'rows','complete');
[corr18(249,1),corr18(249,2)] = corr(PSStotal(:),commRadius(:,3),'rows','complete');
[corr18(250,1),corr18(250,2)] = corr(PSStotal(:),commRadius(:,4),'rows','complete');

%% correlations between PHQ15 and network metrics
[corr18(81,1),corr18(81,2)] = corr(PHQ15total(:),numCommunities(:,2),'rows','complete');
[corr18(82,1),corr18(82,2)] = corr(PHQ15total(:),numCommunities(:,3),'rows','complete');
[corr18(83,1),corr18(83,2)] = corr(PHQ15total(:),numCommunities(:,4),'rows','complete');
[corr18(177,1),corr18(177,2)] = corr(PHQ15total(:),modularity(:,2),'rows','complete');
[corr18(178,1),corr18(178,2)] = corr(PHQ15total(:),modularity(:,3),'rows','complete');
[corr18(179,1),corr18(179,2)] = corr(PHQ15total(:),modularity(:,4),'rows','complete');
[corr18(180,1),corr18(180,2)] = corr(PHQ15total(:),modularity_VA(:,2),'rows','complete');
[corr18(181,1),corr18(181,2)] = corr(PHQ15total(:),modularity_VA(:,3),'rows','complete');
[corr18(182,1),corr18(182,2)] = corr(PHQ15total(:),modularity_VA(:,4),'rows','complete');
[corr18(84,1),corr18(84,2)] = corr(PHQ15total(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(85,1),corr18(85,2)] = corr(PHQ15total(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(86,1),corr18(86,2)] = corr(PHQ15total(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(87,1),corr18(87,2)] = corr(PHQ15total(:),networkDensity(:,3),'rows','complete');
[corr18(88,1),corr18(88,2)] = corr(PHQ15total(:),networkDensity(:,4),'rows','complete');
[corr18(183,1),corr18(183,2)] = corr(PHQ15total(:),networkDensity(:,5),'rows','complete');
[corr18(184,1),corr18(184,2)] = corr(PHQ15total(:),networkDensity(:,6),'rows','complete');
[corr18(185,1),corr18(185,2)] = corr(PHQ15total(:),networkDensity(:,7),'rows','complete');
[corr18(89,1),corr18(89,2)] = corr(PHQ15total(:),participation(:,2),'rows','complete');
[corr18(90,1),corr18(90,2)] = corr(PHQ15total(:),participation(:,3),'rows','complete');
[corr18(91,1),corr18(91,2)] = corr(PHQ15total(:),participation(:,5),'rows','complete');
[corr18(92,1),corr18(92,2)] = corr(PHQ15total(:),participation(:,6),'rows','complete');
[corr18(93,1),corr18(93,2)] = corr(PHQ15total(:),gateway(:,2),'rows','complete'); % node strength
[corr18(94,1),corr18(94,2)] = corr(PHQ15total(:),gateway(:,3),'rows','complete'); % node strength
[corr18(95,1),corr18(95,2)] = corr(PHQ15total(:),gateway(:,5),'rows','complete'); % node strength
[corr18(96,1),corr18(96,2)] = corr(PHQ15total(:),gateway(:,6),'rows','complete'); % node strength
[corr18(97,1),corr18(97,2)] = corr(PHQ15total(:),diversity(:,2),'rows','complete');
[corr18(98,1),corr18(98,2)] = corr(PHQ15total(:),diversity(:,3),'rows','complete');
[corr18(99,1),corr18(99,2)] = corr(PHQ15total(:),diversity(:,5),'rows','complete');
[corr18(100,1),corr18(100,2)] = corr(PHQ15total(:),diversity(:,6),'rows','complete');
[corr18(251,1),corr18(251,2)] = corr(PHQ15total(:),commRadius(:,2),'rows','complete');
[corr18(252,1),corr18(252,2)] = corr(PHQ15total(:),commRadius(:,3),'rows','complete');
[corr18(253,1),corr18(253,2)] = corr(PHQ15total(:),commRadius(:,4),'rows','complete');

%% correlations between PHQ8 and network metrics
[corr18(101,1),corr18(101,2)] = corr(PHQ8total(:),numCommunities(:,2),'rows','complete');
[corr18(102,1),corr18(102,2)] = corr(PHQ8total(:),numCommunities(:,3),'rows','complete');
[corr18(103,1),corr18(103,2)] = corr(PHQ8total(:),numCommunities(:,4),'rows','complete');
[corr18(186,1),corr18(186,2)] = corr(PHQ8total(:),modularity(:,2),'rows','complete');
[corr18(187,1),corr18(187,2)] = corr(PHQ8total(:),modularity(:,3),'rows','complete');
[corr18(188,1),corr18(188,2)] = corr(PHQ8total(:),modularity(:,4),'rows','complete');
[corr18(189,1),corr18(189,2)] = corr(PHQ8total(:),modularity_VA(:,2),'rows','complete');
[corr18(190,1),corr18(190,2)] = corr(PHQ8total(:),modularity_VA(:,3),'rows','complete');
[corr18(191,1),corr18(191,2)] = corr(PHQ8total(:),modularity_VA(:,4),'rows','complete');
[corr18(104,1),corr18(104,2)] = corr(PHQ8total(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(105,1),corr18(105,2)] = corr(PHQ8total(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(106,1),corr18(106,2)] = corr(PHQ8total(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(107,1),corr18(107,2)] = corr(PHQ8total(:),networkDensity(:,3),'rows','complete');
[corr18(108,1),corr18(108,2)] = corr(PHQ8total(:),networkDensity(:,4),'rows','complete');
[corr18(192,1),corr18(192,2)] = corr(PHQ8total(:),networkDensity(:,5),'rows','complete');
[corr18(193,1),corr18(193,2)] = corr(PHQ8total(:),networkDensity(:,6),'rows','complete');
[corr18(194,1),corr18(194,2)] = corr(PHQ8total(:),networkDensity(:,7),'rows','complete');
[corr18(109,1),corr18(109,2)] = corr(PHQ8total(:),participation(:,2),'rows','complete');
[corr18(110,1),corr18(110,2)] = corr(PHQ8total(:),participation(:,3),'rows','complete');
[corr18(111,1),corr18(111,2)] = corr(PHQ8total(:),participation(:,5),'rows','complete');
[corr18(112,1),corr18(112,2)] = corr(PHQ8total(:),participation(:,6),'rows','complete');
[corr18(113,1),corr18(113,2)] = corr(PHQ8total(:),gateway(:,2),'rows','complete'); % node strength
[corr18(114,1),corr18(114,2)] = corr(PHQ8total(:),gateway(:,3),'rows','complete'); % node strength
[corr18(115,1),corr18(115,2)] = corr(PHQ8total(:),gateway(:,5),'rows','complete'); % node strength
[corr18(116,1),corr18(116,2)] = corr(PHQ8total(:),gateway(:,6),'rows','complete'); % node strength
[corr18(117,1),corr18(117,2)] = corr(PHQ8total(:),diversity(:,2),'rows','complete');
[corr18(118,1),corr18(118,2)] = corr(PHQ8total(:),diversity(:,3),'rows','complete');
[corr18(119,1),corr18(119,2)] = corr(PHQ8total(:),diversity(:,5),'rows','complete');
[corr18(120,1),corr18(120,2)] = corr(PHQ8total(:),diversity(:,6),'rows','complete');
[corr18(254,1),corr18(254,2)] = corr(PHQ8total(:),commRadius(:,2),'rows','complete');
[corr18(255,1),corr18(255,2)] = corr(PHQ8total(:),commRadius(:,3),'rows','complete');
[corr18(256,1),corr18(256,2)] = corr(PHQ8total(:),commRadius(:,4),'rows','complete');

%% correlations between neuroticism and network metrics
[corr18(121,1),corr18(121,2)] = corr(neuroticism(:),numCommunities(:,2),'rows','complete');
[corr18(122,1),corr18(122,2)] = corr(neuroticism(:),numCommunities(:,3),'rows','complete');
[corr18(123,1),corr18(123,2)] = corr(neuroticism(:),numCommunities(:,4),'rows','complete');
[corr18(195,1),corr18(195,2)] = corr(neuroticism(:),modularity(:,2),'rows','complete');
[corr18(196,1),corr18(196,2)] = corr(neuroticism(:),modularity(:,3),'rows','complete');
[corr18(197,1),corr18(197,2)] = corr(neuroticism(:),modularity(:,4),'rows','complete');
[corr18(198,1),corr18(198,2)] = corr(neuroticism(:),modularity_VA(:,2),'rows','complete');
[corr18(199,1),corr18(199,2)] = corr(neuroticism(:),modularity_VA(:,3),'rows','complete');
[corr18(200,1),corr18(200,2)] = corr(neuroticism(:),modularity_VA(:,4),'rows','complete');
[corr18(124,1),corr18(124,2)] = corr(neuroticism(:),clusterCoefGlobal(:,2),'rows','complete');
[corr18(125,1),corr18(125,2)] = corr(neuroticism(:),clusterCoefGlobal(:,3),'rows','complete');
[corr18(126,1),corr18(126,2)] = corr(neuroticism(:),clusterCoefGlobal(:,4),'rows','complete');
[corr18(127,1),corr18(127,2)] = corr(neuroticism(:),networkDensity(:,3),'rows','complete');
[corr18(128,1),corr18(128,2)] = corr(neuroticism(:),networkDensity(:,4),'rows','complete');
[corr18(201,1),corr18(201,2)] = corr(neuroticism(:),networkDensity(:,5),'rows','complete');
[corr18(202,1),corr18(202,2)] = corr(neuroticism(:),networkDensity(:,6),'rows','complete');
[corr18(203,1),corr18(203,2)] = corr(neuroticism(:),networkDensity(:,7),'rows','complete');
[corr18(129,1),corr18(129,2)] = corr(neuroticism(:),participation(:,2),'rows','complete');
[corr18(130,1),corr18(130,2)] = corr(neuroticism(:),participation(:,3),'rows','complete');
[corr18(131,1),corr18(131,2)] = corr(neuroticism(:),participation(:,5),'rows','complete');
[corr18(132,1),corr18(132,2)] = corr(neuroticism(:),participation(:,6),'rows','complete');
[corr18(133,1),corr18(133,2)] = corr(neuroticism(:),gateway(:,2),'rows','complete'); % node strength
[corr18(134,1),corr18(134,2)] = corr(neuroticism(:),gateway(:,3),'rows','complete'); % node strength
[corr18(135,1),corr18(135,2)] = corr(neuroticism(:),gateway(:,5),'rows','complete'); % node strength
[corr18(136,1),corr18(136,2)] = corr(neuroticism(:),gateway(:,6),'rows','complete'); % node strength
[corr18(137,1),corr18(137,2)] = corr(neuroticism(:),diversity(:,2),'rows','complete');
[corr18(138,1),corr18(138,2)] = corr(neuroticism(:),diversity(:,3),'rows','complete');
[corr18(139,1),corr18(139,2)] = corr(neuroticism(:),diversity(:,5),'rows','complete');
[corr18(140,1),corr18(140,2)] = corr(neuroticism(:),diversity(:,6),'rows','complete');
[corr18(257,1),corr18(257,2)] = corr(neuroticism(:),commRadius(:,2),'rows','complete');
[corr18(258,1),corr18(258,2)] = corr(neuroticism(:),commRadius(:,3),'rows','complete');
[corr18(259,1),corr18(259,2)] = corr(neuroticism(:),commRadius(:,4),'rows','complete');

%% correlations between zICC and psychological measures
[corr18(204,1),corr18(204,2)] = corr(ASItotal(:),rtoZ(:,3),'rows','complete');
[corr18(205,1),corr18(205,2)] = corr(ERStotal(:),rtoZ(:,3),'rows','complete');
[corr18(206,1),corr18(206,2)] = corr(GADtotal(:),rtoZ(:,3),'rows','complete');
[corr18(207,1),corr18(207,2)] = corr(PSStotal(:),rtoZ(:,3),'rows','complete');
[corr18(208,1),corr18(208,2)] = corr(PHQ15total(:),rtoZ(:,3),'rows','complete');
[corr18(209,1),corr18(209,2)] = corr(PHQ8total(:),rtoZ(:,3),'rows','complete');
[corr18(210,1),corr18(210,2)] = corr(neuroticism(:),rtoZ(:,3),'rows','complete');

%% correlations between affect and psychological measures
[corr18(211,1),corr18(211,2)] = corr(ASItotal(:),affect(:,1),'rows','complete');
[corr18(212,1),corr18(212,2)] = corr(ASItotal(:),affect(:,2),'rows','complete');
[corr18(213,1),corr18(213,2)] = corr(ASItotal(:),affect(:,3),'rows','complete');
[corr18(214,1),corr18(214,2)] = corr(ASItotal(:),affect(:,4),'rows','complete');
[corr18(215,1),corr18(215,2)] = corr(ERStotal(:),affect(:,1),'rows','complete');
[corr18(216,1),corr18(216,2)] = corr(ERStotal(:),affect(:,2),'rows','complete');
[corr18(217,1),corr18(217,2)] = corr(ERStotal(:),affect(:,3),'rows','complete');
[corr18(218,1),corr18(218,2)] = corr(ERStotal(:),affect(:,4),'rows','complete');
[corr18(219,1),corr18(219,2)] = corr(GADtotal(:),affect(:,1),'rows','complete');
[corr18(220,1),corr18(220,2)] = corr(GADtotal(:),affect(:,2),'rows','complete');
[corr18(221,1),corr18(221,2)] = corr(GADtotal(:),affect(:,3),'rows','complete');
[corr18(222,1),corr18(222,2)] = corr(GADtotal(:),affect(:,4),'rows','complete');
[corr18(223,1),corr18(223,2)] = corr(PSStotal(:),affect(:,1),'rows','complete');
[corr18(224,1),corr18(224,2)] = corr(PSStotal(:),affect(:,2),'rows','complete');
[corr18(225,1),corr18(225,2)] = corr(PSStotal(:),affect(:,3),'rows','complete');
[corr18(226,1),corr18(226,2)] = corr(PSStotal(:),affect(:,4),'rows','complete');
[corr18(227,1),corr18(227,2)] = corr(PHQ15total(:),affect(:,1),'rows','complete');
[corr18(228,1),corr18(228,2)] = corr(PHQ15total(:),affect(:,2),'rows','complete');
[corr18(229,1),corr18(229,2)] = corr(PHQ15total(:),affect(:,3),'rows','complete');
[corr18(230,1),corr18(230,2)] = corr(PHQ15total(:),affect(:,4),'rows','complete');
[corr18(231,1),corr18(231,2)] = corr(PHQ8total(:),affect(:,1),'rows','complete');
[corr18(232,1),corr18(232,2)] = corr(PHQ8total(:),affect(:,2),'rows','complete');
[corr18(233,1),corr18(233,2)] = corr(PHQ8total(:),affect(:,3),'rows','complete');
[corr18(234,1),corr18(234,2)] = corr(PHQ8total(:),affect(:,4),'rows','complete');
[corr18(235,1),corr18(235,2)] = corr(neuroticism(:),affect(:,1),'rows','complete');
[corr18(236,1),corr18(236,2)] = corr(neuroticism(:),affect(:,2),'rows','complete');
[corr18(237,1),corr18(237,2)] = corr(neuroticism(:),affect(:,3),'rows','complete');
[corr18(238,1),corr18(238,2)] = corr(neuroticism(:),affect(:,4),'rows','complete');

%% organize and save output
correlations = array2table(corr18);
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'v1','v2','gamma','threshold','r','p'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_criterion_validity_checks_polychoric.xlsx'];
    writetable(correlations,filename);
end