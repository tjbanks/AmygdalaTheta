# Upload analysis
import slack
import os

def upload(text='message text',filename='analysis.png', channel='general'):
    token = os.getenv("SLACK_TOKEN") # yay security
    if not token:
        print("set your SLACK_TOKEN environment variable first")
        return
    client = slack.WebClient(token=token)
    img = open(filename, 'rb').read()
    new_file = client.files_upload(channels='C051KAS1PHN', initial_comment='', filename='analysis.png', content=img)
    image_file_path = new_file.get('file').get('permalink')
    print(image_file_path)
    attachments = [{"title": "analysis.png", "image_url": image_file_path}]
    client.chat_postMessage(channel=channel, text=text)

if __name__ == '__main__':
    upload()
