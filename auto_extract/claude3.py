import anthropic
from openai import OpenAI
from pypdf import PdfReader
from tqdm import tqdm
from abc import ABC, abstractmethod

from .config import Config

class BaseClient(ABC):
    def __init__(self, pdf_path: str = None):
        self.api_key = Config().API_KEY
        self.api_base_url = Config().API_BASE_URL
        self.model = Config().MODEL
        self.messages = []
        
        self.pdf_path = pdf_path
        self.extracted_text = ''
        if self.pdf_path:
            if self.pdf_path.endswith('.txt'):
                with open(self.pdf_path, 'r') as f:
                    self.extracted_text = f.read()
            else:
                self.extracted_text = self.extract_text_from_pdf(self.pdf_path)
    
    def extract_text_from_pdf(self, pdf_path: str) -> str:
        reader = PdfReader(pdf_path)
        text = ''
        for page in tqdm(reader.pages, desc='Extracting text from PDF'):
            text += page.extract_text()
            text += '\n'
        return text
    
    def update_messages(func):
        def wrapper(self, *args, **kwargs):
            resp = func(self, *args, **kwargs)
            self.messages.append({"role": "assistant", "content": resp})
            return resp
        return wrapper
    
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
    def chat(self, new_message: str, max_tokens = 1000):

        self.messages.append({"role": "user", "content": new_message})
        resp = self._retrieve(self.messages, max_tokens=max_tokens)
        
        return resp
    
    @abstractmethod
    def create_client(self):
        pass
    
    @abstractmethod
    def _retrieve(self):
        pass

class Openai(BaseClient):
    def __init__(self, pdf_path: str = None):
        super().__init__(pdf_path)
        self.create_client()

    def create_client(self):
        self.client = OpenAI(api_key=self.api_key, api_base_url=self.api_base_url)

    # add retrieve method   
    def _retrieve(self, 
                 messages: list, 
                 max_tokens: int = 1000):
        
        completion = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens
        )
        
        return completion.choices[0].message.content
    
class Claude3(BaseClient):
    def __init__(self, pdf_path: str = None):
        super().__init__(pdf_path)
        self.create_client()
        
    def create_client(self):
        self.client = anthropic.Anthropic(api_key=self.api_key)
    
    # add retrieve method
    def _retrieve(self, 
                 messages: list, 
                 max_tokens: int = 1000):
        
        completion = self.client.messages.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens
        )
        
        return completion.content[0].text

    def update_messages(func):
        def wrapper(self, *args, **kwargs):
            resp = func(self, *args, **kwargs)
            self.messages.append({"role": "assistant", "content": resp})
            return resp
        return wrapper

    @update_messages
    def initiate_propmt(self, 
                        initial_message: str = Config().get_prompt('SYSTEM_PROMPT'),
                        user_message: str = Config().get_prompt('USER_ARTICLE_PROMPT'),
                        ):
        
        if self.extracted_text == '':
            raise ValueError('No PDF file provided')

         # claude3 does not have a system message
        initial_message += user_message
        initial_message += self.extracted_text
        self.messages.append({"role": "user", "content": initial_message})
        
        resp = self._retrieve(self.messages)
        
        return resp