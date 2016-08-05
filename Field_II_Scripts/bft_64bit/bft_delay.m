%BFT_DELAY Apply a delay on one line.
%   This function is not particularly useful, unless someones wants
%   to make phase-aberration or motion compensation. Even in those 
%   cases, this is not a "good" function to use, since the changes
%   are not applied during a beamformation procedure. 
%   The function however is capable of applying different delays
%   in the different image zones, and therefore is more powerful and
%   faster than a simple shift of an array in Matlab.
%
%USAGE : [out_line] = bft_delay(in_line, delays, times, input_start_time...
%                          output_start_time, no_output_samples, use_filter)
%
%INPUTS: in_line - A vector with the samples
%        delays  - A vector with the delays to be applied. 
%        times   - A vector saying when the next delay is valid. The first
%                  delay is applied until times(2). delays(2) is applied 
%                  after times(2), etc.
%        input_start_time - The starting time of the first sample in the 
%                  input signal.
%        output_start_time - When is the first output sample to start.
%        no_output_samples - Number of samples in "out_line". 
%        use_filter - Whether to use a linear phase filter or not.
%                     This one is optional. The default is not to use it.
%                    
%
%OUTPUT: out_line - The delayed line
%
%CREATED: 07 Sep 2000, Svetoslav Nikolov
%

function out_line = bft_delay(in_line, delays, times, input_start_time, ...
                              output_start_time, no_output_samples, ...
                              use_filter)

if (isempty(use_filter)) use_filter = 0; end;

[m n] = size(in_line);
if (n ==1)
   if (use_filter == 0)
       out_line = bft(19, in_line, delays, times, input_start_time, ...
                                  output_start_time, no_output_samples);
   else
       out_line = bft(20, in_line, delays, times, input_start_time, ...
                                  output_start_time, no_output_samples);

   end                              
else
   [k, l] = size(delays)
   if (k == 1)
      if (use_filter == 0)
         for ii = 1:n
            out_line(1:no_output_samples,ii) =  bft(19, in_line(:,ii), delays, times, input_start_time, ...
                               output_start_time, no_output_samples);
         end
      else
         for ii = 1:n
             out_line(1:no_output_samples,ii) = bft(20, in_line(:,ii), delays, times, input_start_time, ...
                                  output_start_time, no_output_samples);
         end
      end
   else
      [g h] = size(times);
      if (g ~= k | h~=l)
         error('The size of times and delays must be the same');
      end
      if (l~=n)
         error('You need one column if delays per one column of data ');
      end
      if (use_filter == 0)
         for ii = 1:n
            out_line(1:no_output_samples,ii) =  bft(19, in_line(:,ii), delays(:,ii), times(:,ii), input_start_time, ...
                               output_start_time, no_output_samples);
         end
      else
         for ii = 1:n
				 
             out_line(1:no_output_samples,ii) = bft(20, in_line(:,ii), delays(:,ii), times(:,ii), input_start_time, ...
                                  output_start_time, no_output_samples);
         end
      end
   end
end
