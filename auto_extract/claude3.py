# import requests
from openai import OpenAI
from pypdf import PdfReader
from tqdm import tqdm

from .config import Config

class Claude3:
    def __init__(self, pdf_path: str = None):
        self.api_key = Config().API_KEY
        self.api_base_url = Config().API_BASE_URL
        self.client = OpenAI(api_key=self.api_key, base_url=self.api_base_url)
        self.messages = []
        
        self.pdf_path = pdf_path
        self.extracted_text = ''
        if self.pdf_path:
            if self.pdf_path.endswith('.txt'):
                with open(self.pdf_path, 'r') as f:
                    self.extracted_text = f.read()
            else:
                self.extracted_text = self.extract_text_from_pdf(self.pdf_path)
            
    def update_messages(func):
        def wrapper(self, *args, **kwargs):
            resp = func(self, *args, **kwargs)
            self.messages.append({"role": "assistant", "content": resp})
            return resp
        return wrapper

    def extract_text_from_pdf(self, pdf_path: str) -> str:
        reader = PdfReader(pdf_path)
        text = ''
        for page in tqdm(reader.pages, desc='Extracting text from PDF'):
            text += page.extract_text()
            text += '\n'
        return text
    
    # add retrieve method   
    def _retrieve(self, 
                 messages: list, 
                 max_tokens: int = 1000):
        
        completion = self.client.chat.completions.create(
            model="claude-3-sonnet-20240229",
            messages=messages,
            max_tokens=max_tokens
        )
        
        return completion.choices[0].message.content
    
    @update_messages
    def initiate_propmt(self, 
                        initial_message: str = Config().get_prompt('SYSTEM_PROMPT'),
                        user_message: str = Config().get_prompt('USER_ARTICLE_PROMPT'),
                        ):
        
        if self.extracted_text == '':
            raise ValueError('No PDF file provided')

        self.messages.append({"role": "system", "content": initial_message})
        user_message += self.extracted_text
        self.messages.append({"role": "user", "content": user_message})
        
        resp = self._retrieve(self.messages)
        
        return resp
    
    @update_messages
    def chat(self, new_message: str):

        self.messages.append({"role": "user", "content": new_message})
        resp = self._retrieve(self.messages)
        
        return resp
    
