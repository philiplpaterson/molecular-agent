import { useState } from "react";
import ChatWindow from "./components/Chat/ChatWindow";

function App() {
  const [conversationId, setConversationId] = useState<string | null>(null);

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
      <header className="bg-white shadow-sm border-b border-gray-200">
        <div className="max-w-4xl mx-auto px-4 py-4">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 bg-gradient-to-br from-blue-500 to-indigo-600 rounded-lg flex items-center justify-center">
              <span className="text-white text-xl">🧬</span>
            </div>
            <div>
              <h1 className="text-xl font-bold text-gray-900">MolecularAgent</h1>
              <p className="text-sm text-gray-500">AI Drug Discovery Assistant</p>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-4xl mx-auto px-4 py-6">
        <ChatWindow
          conversationId={conversationId}
          onConversationIdChange={setConversationId}
        />
      </main>

      <footer className="fixed bottom-0 left-0 right-0 bg-white border-t border-gray-200 py-2">
        <div className="max-w-4xl mx-auto px-4 text-center text-xs text-gray-500">
          Powered by PydanticAI and RDKit
        </div>
      </footer>
    </div>
  );
}

export default App;
