%% Send Email after finishing
function send_mail_message(id, mail, subject, password, message,attachment)
%% SEND_MAIL_MESSAGE send email to gmail once calculation is done
% Example
% send_mail_message('its.neeraj','Simulation finished','This is the message area','results.doc')
 
% Pradyumna
% June 2008

% Modified by Julian Hinz, @Luthi lab 2018

%% Set up Gmail SMTP service.
% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% Gmail server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% Send the email

sendmail(id,subject,message,attachment)

end