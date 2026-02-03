import ReactMarkdown from "react-markdown";
import ToolCallDisplay from "./ToolCallDisplay";

export interface ToolCall {
  name: string;
  arguments: Record<string, unknown>;
  result: Record<string, unknown>;
}

export interface Message {
  id: string;
  role: "user" | "assistant";
  content: string;
  toolCalls?: ToolCall[];
}

interface MessageListProps {
  messages: Message[];
}

export default function MessageList({ messages }: MessageListProps) {
  return (
    <div className="space-y-4">
      {messages.map((message) => (
        <div
          key={message.id}
          className={`flex ${
            message.role === "user" ? "justify-end" : "justify-start"
          }`}
        >
          <div
            className={`max-w-[85%] ${
              message.role === "user"
                ? "bg-blue-600 text-white rounded-2xl rounded-br-md px-4 py-2"
                : "bg-gray-100 text-gray-900 rounded-2xl rounded-bl-md px-4 py-3"
            }`}
          >
            {message.role === "assistant" ? (
              <div className="space-y-3">
                {message.toolCalls && message.toolCalls.length > 0 && (
                  <div className="space-y-2">
                    {message.toolCalls.map((toolCall, index) => (
                      <ToolCallDisplay key={index} toolCall={toolCall} />
                    ))}
                  </div>
                )}
                <div className="prose prose-sm max-w-none">
                  <ReactMarkdown>{message.content}</ReactMarkdown>
                </div>
              </div>
            ) : (
              <p>{message.content}</p>
            )}
          </div>
        </div>
      ))}
    </div>
  );
}
