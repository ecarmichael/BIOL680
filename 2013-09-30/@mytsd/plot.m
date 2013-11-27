function [h] = plot(tsd_in,varargin)
% function plot(tsd_in)
%
% plots the data in a tsd object
tsd_in.t = Range(tsd_in);
tsd_in.data = Data(tsd_in);
 
h = plot(tsd_in.t,tsd_in.data,varargin{:});
